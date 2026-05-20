"""Writer for the NetCDF Rupture Format (NRF).

The NRF is the binary form of a kinematic rupture model consumed by SeisSol.
This module produces files that match the layout written by the reference C++
``rconv`` tool closely enough that SeisSol's NRF reader treats them as
equivalent. See ``preprocessing/science/rconv/spec/nrf.cdl`` in the SeisSol
repository for the canonical schema.

A note on units
---------------

SRF carries lengths in cm and times in s; the NRF schema requires SI lengths
(metres). This module performs the conversions explicitly:

* depth (km) → m  (via :class:`rconv.mapping.CoordinateMapping`)
* area  (cm²) → m²       (× 1e-4)
* shear modulus (g/(cm·s²)) → Pa = kg/(m·s²)  (× 1e-1)
* slip rate (cm/s) → m/s (× 1e-2)

A note on the ``Subfault_units`` compound
-----------------------------------------

The C++ ``rconv`` attaches a compound-typed ``units`` attribute on
``subfaults`` (compound of seven VLEN strings). The Python ``netCDF4`` package
does not support strings inside compound types, so we attach the same
information as a plain-text attribute instead. SeisSol's reader never inspects
this attribute — it's pure metadata for tools like ``ncdump``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Sequence

import netCDF4
import numpy as np

from .mapping import CoordinateMapping
from .srf import PointSource

_log = logging.getLogger(__name__)

# Unit-conversion constants (cgs → SI).
_CM2_TO_M2 = 1.0e-4
_CGS_MU_TO_PA = 1.0e-1
_CMS_TO_MS = 1.0e-2

# numpy dtypes matching the NRF compound types.
_VECTOR3_DTYPE = np.dtype([("x", "f8"), ("y", "f8"), ("z", "f8")])
_SUBFAULT_DTYPE = np.dtype([
    ("tinit",    "f8"),
    ("timestep", "f8"),
    ("mu",       "f8"),
    ("area",     "f8"),
    ("tan1",     _VECTOR3_DTYPE),
    ("tan2",     _VECTOR3_DTYPE),
    ("normal",   _VECTOR3_DTYPE),
])

# Units metadata, mirroring rconv's Subfault_units compound attribute.
_SUBFAULT_UNITS_TEXT = (
    "tinit=s, timestep=s, mu=pascal, area=m^2, tan1=m, tan2=m, normal=m"
)


def write_nrf(
    path: str | Path,
    sources: Sequence[PointSource],
    mapping: CoordinateMapping,
    *,
    normalize_onset: bool = False,
) -> None:
    """Write a list of point sources to a NetCDF Rupture Format file.

    Parameters
    ----------
    path:
        Destination file path. Overwritten if it exists.
    sources:
        Point sources, as produced by :func:`rconv.srf.parse_srf`. Must be
        non-empty.
    mapping:
        Coordinate mapping into the mesh CS.
    normalize_onset:
        If true, subtract the minimum ``tinit`` across all sources from every
        source's onset. Useful when the SRF's absolute clock is meaningless.
    """
    if not sources:
        raise ValueError("write_nrf: empty source list")

    centres, subfaults, sroffsets, sliprates = _build_arrays(
        sources, mapping, normalize_onset=normalize_onset,
    )

    path = Path(path)
    # mode='w' clobbers — matches NC_CLOBBER in the C++ writer.
    with netCDF4.Dataset(path, mode="w", format="NETCDF4") as ds:
        _write_dataset(ds, centres, subfaults, sroffsets, sliprates)


# ---------------------------------------------------------------------------
# array assembly
# ---------------------------------------------------------------------------

def _build_arrays(
    sources: Sequence[PointSource],
    mapping: CoordinateMapping,
    *,
    normalize_onset: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    n = len(sources)

    # Offset adjustment for onsets.
    if normalize_onset:
        min_tinit = min(s.tinit for s in sources)
        _log.info("Normalising onset: minimal tinit %g s → 0.", min_tinit)
    else:
        min_tinit = 0.0

    # Allocate output arrays.
    centres = np.empty(n, dtype=_VECTOR3_DTYPE)
    subfaults = np.empty(n, dtype=_SUBFAULT_DTYPE)

    # sroffsets layout: (n+1, 3). offsets[i][d] is the start index of source i
    # in direction d's flat sliprate array; offsets[n][d] is the total length.
    sroffsets = np.zeros((n + 1, 3), dtype=np.uint32)

    # Accumulate sliprate samples per direction.
    sliprate_chunks: tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]] = ([], [], [])

    for i, src in enumerate(sources):
        x, y, z = mapping.to_mesh(src.longitude, src.latitude, src.depth).as_tuple()
        centres[i] = (x, y, z)

        # Three orthogonal unit vectors of the fault-local frame, transformed
        # into the MCS. These define the slip basis in the mesh CS.
        tan1 = mapping.to_mcs(src.strike, src.dip, src.rake, 1.0, 0.0, 0.0).as_tuple()
        tan2 = mapping.to_mcs(src.strike, src.dip, src.rake, 0.0, 1.0, 0.0).as_tuple()
        normal = mapping.to_mcs(src.strike, src.dip, src.rake, 0.0, 0.0, 1.0).as_tuple()

        subfaults[i] = (
            src.tinit - min_tinit,                  # tinit (s)
            src.dt,                                 # timestep (s)
            src.shear_modulus * _CGS_MU_TO_PA,      # mu (Pa)
            src.area * _CM2_TO_M2,                  # area (m^2)
            tan1, tan2, normal,
        )

        # Sliprate offsets: cumulative count per direction.
        for d in range(3):
            samples = np.asarray(src.slip_rates[d], dtype=np.float64) * _CMS_TO_MS
            sliprate_chunks[d].append(samples)
            sroffsets[i + 1, d] = sroffsets[i, d] + samples.size

    # Concatenate per direction. Empty list → zero-length array.
    sliprates: tuple[np.ndarray, np.ndarray, np.ndarray] = tuple(
        np.concatenate(chunks) if chunks else np.empty(0, dtype=np.float64)
        for chunks in sliprate_chunks
    )  # type: ignore[assignment]

    return centres, subfaults, sroffsets, sliprates


# ---------------------------------------------------------------------------
# NetCDF emission
# ---------------------------------------------------------------------------

def _write_dataset(
    ds: netCDF4.Dataset,
    centres: np.ndarray,
    subfaults: np.ndarray,
    sroffsets: np.ndarray,
    sliprates: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    # --- compound types ---------------------------------------------------
    vec3_t = ds.createCompoundType(_VECTOR3_DTYPE, "Vector3")
    sf_t = ds.createCompoundType(_SUBFAULT_DTYPE, "Subfault")

    # --- dimensions -------------------------------------------------------
    n = centres.shape[0]
    ds.createDimension("source", n)
    ds.createDimension("sroffset", n + 1)
    ds.createDimension("direction", 3)
    # netCDF4-python maps size=0 to UNLIMITED; this is a cosmetic difference
    # from the C++ writer (which makes them fixed-size-0). The reader doesn't
    # care.
    ds.createDimension("sample1", sliprates[0].size)
    ds.createDimension("sample2", sliprates[1].size)
    ds.createDimension("sample3", sliprates[2].size)

    # --- variables --------------------------------------------------------
    centres_var = ds.createVariable("centres", vec3_t, ("source",))
    centres_var.units = "m"
    centres_var[:] = centres

    subfaults_var = ds.createVariable(
        "subfaults", sf_t, ("source",), zlib=True, complevel=1, shuffle=False,
    )
    subfaults_var.units = _SUBFAULT_UNITS_TEXT
    subfaults_var[:] = subfaults

    sroffsets_var = ds.createVariable(
        "sroffsets", "u4", ("sroffset", "direction"),
        zlib=True, complevel=1, shuffle=False,
    )
    sroffsets_var[:] = sroffsets

    for d, name in enumerate(("sliprates1", "sliprates2", "sliprates3")):
        dim = f"sample{d + 1}"
        var = ds.createVariable(
            name, "f8", (dim,), zlib=True, complevel=1, shuffle=False,
        )
        var.units = "m/s"
        if sliprates[d].size:
            var[:] = sliprates[d]


# ---------------------------------------------------------------------------
# reader (for tests + round-trip verification)
# ---------------------------------------------------------------------------

def read_nrf(path: str | Path) -> dict:
    """Read an NRF file and return a dict of its arrays.

    Provided primarily for testing and round-trip verification; not intended
    as a replacement for SeisSol's own consumer. The return value keys mirror
    the variable names in the NRF schema: ``centres``, ``subfaults``,
    ``sroffsets``, ``sliprates1/2/3``.
    """
    path = Path(path)
    with netCDF4.Dataset(path, mode="r") as ds:
        return {
            "centres":    np.asarray(ds["centres"][:]),
            "subfaults":  np.asarray(ds["subfaults"][:]),
            "sroffsets":  np.asarray(ds["sroffsets"][:]),
            "sliprates1": np.asarray(ds["sliprates1"][:]),
            "sliprates2": np.asarray(ds["sliprates2"][:]),
            "sliprates3": np.asarray(ds["sliprates3"][:]),
        }
