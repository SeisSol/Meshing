"""rconv-py — convert SRF kinematic rupture models to SeisSol's NRF.

A pure-Python replacement for the C++ ``rconv`` tool that ships with SeisSol.
This package is dependency-free with respect to SeisSol; it relies only on
``numpy``, ``pyproj`` and ``netCDF4``.

The public API is small:

* :func:`rconv.srf.parse_srf` — read an SRF file into :class:`PointSource`\\ s.
* :class:`rconv.mapping.CoordinateMapping` — transform geographic coords and
  fault-local vectors into a mesh CS.
* :func:`rconv.nrf.write_nrf` — emit a NetCDF Rupture Format file for SeisSol.
* :func:`rconv.xmf.write_xmf` — emit an XDMF file for ParaView.

Example
-------

>>> from rconv import parse_srf, CoordinateMapping, write_nrf
>>> sources = parse_srf("northridge.srf")
>>> mapping = CoordinateMapping(
...     "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=-118.515 +lat_0=34.344 "
...     "+units=m +axis=ned"
... )
>>> write_nrf("northridge.nrf", sources, mapping)
"""

from __future__ import annotations

from .mapping import CoordinateMapping, Vec3
from .nrf import read_nrf, write_nrf
from .srf import PointSource, SRFParseError, parse_srf
from .xmf import write_xmf

__all__ = [
    "CoordinateMapping",
    "PointSource",
    "SRFParseError",
    "Vec3",
    "parse_srf",
    "read_nrf",
    "write_nrf",
    "write_xmf",
]

__version__ = "0.1.0"
