"""Coordinate transformations for SRF → NRF conversion.

Two distinct operations live here:

1. ``to_mesh`` — transform a geographic point (lon, lat, depth) into the user's
   mesh coordinate system (MCS) using a PROJ.4 / PROJ string.
2. ``to_mcs`` — rotate a vector given in the fault-local frame
   (along strike, along dip, along normal) into the MCS.

Both honour the ``+axis=`` PROJ option (default ``enu``: x=east, y=north, z=up).
The axis option determines how a canonical (north, east, down) triple — which is
the natural frame of the SRF format — is permuted/signed into the user's MCS.

A note on pyproj subtleties
---------------------------

pyproj's behaviour around ``+axis=`` differs between projected and geocentric
target CRSs:

* For **geocentric** targets, pyproj fully honours ``+axis=`` (it permutes and
  negates X/Y/Z accordingly), so we delegate the 3-D transform to it.
* For **projected** targets, pyproj only permutes the *horizontal* x/y when
  ``always_xy=False``; the vertical z is always passed through. To get
  consistent, predictable behaviour we therefore do the horizontal projection
  with ``always_xy=True`` (yielding canonical east/north) and apply the axis
  permutation ourselves.

The end result matches the C++ ``rconv``'s use of PROJ.4's ``pj_transform`` plus
the ``adjustAxes`` helper, without depending on PROJ.4's deprecated internals.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass

from pyproj import CRS, Transformer

_log = logging.getLogger(__name__)

# Canonical source CRS for SRF inputs. We do unit conversions (km→m) and
# vertical sign handling ourselves, so the source CRS is the plain WGS84
# geographic CRS in degrees with no vertical hint.
_DEFAULT_SOURCE_PROJ = "+proj=longlat +datum=WGS84 +no_defs"

# Default +axis= in PROJ.4 if the option is absent.
_DEFAULT_AXIS = "enu"

_VALID_AXIS_CHARS = frozenset("enwsud")


@dataclass(frozen=True)
class Vec3:
    """A tiny 3-vector used purely for clarity at call sites."""
    x: float
    y: float
    z: float

    def as_tuple(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


class CoordinateMapping:
    """Transformation between WGS84 lon/lat/depth and a target mesh CS.

    Parameters
    ----------
    target_proj:
        PROJ string describing the MCS, e.g.
        ``"+proj=utm +zone=10 +datum=WGS84 +units=m +axis=ned +no_defs"``.
    source_proj:
        PROJ string for the input CRS. Defaults to plain WGS84 longlat
        (``+proj=longlat +datum=WGS84 +no_defs``); this matches the C++ tool's
        intent, with the km→m vertical conversion handled explicitly in
        :meth:`to_mesh`.
    """

    def __init__(self, target_proj: str, source_proj: str | None = None):
        self.target_proj = target_proj
        self.source_proj = source_proj or _DEFAULT_SOURCE_PROJ

        self.target_crs = CRS.from_proj4(target_proj)
        self.source_crs = CRS.from_proj4(self.source_proj)

        self.axis = _parse_axis(target_proj)
        if len(self.axis) != 3 or any(c not in _VALID_AXIS_CHARS for c in self.axis):
            raise ValueError(
                f"Invalid +axis= in {target_proj!r}: {self.axis!r} "
                f"(expected a 3-character string of {sorted(_VALID_AXIS_CHARS)})"
            )

        self.is_geocentric = self.target_crs.is_geocentric

        # We use always_xy=True for projected targets so we get canonical
        # (east, north) and apply our own axis logic. For geocentric, the
        # value of always_xy is moot — pyproj handles +axis= directly.
        self._transformer = Transformer.from_crs(
            self.source_crs, self.target_crs, always_xy=True,
        )

        # Sanity-check vertical units. The output is documented as meters; if
        # the MCS chose something else, warn (the C++ tool emits a similar
        # warning, since downstream code assumes m / m^2 / m/s).
        self._maybe_warn_units()

    # ------------------------------------------------------------------
    # main operations
    # ------------------------------------------------------------------

    def to_mesh(self, longitude: float, latitude: float, depth_km: float) -> Vec3:
        """Map (lon, lat, depth) → (x, y, z) in the mesh CS.

        ``depth_km`` is positive *downwards*, as in SRF. The returned ``z``
        sign follows the third character of ``self.axis`` (``u``: positive
        upwards → output negative for positive depth; ``d``: positive downwards
        → output positive for positive depth).
        """
        depth_m = depth_km * 1000.0

        if self.is_geocentric:
            # For geocentric targets pyproj does the full 3-D transformation,
            # including any +axis= permutation. We pass height = -depth.
            x, y, z = self._transformer.transform(longitude, latitude, -depth_m)
            return Vec3(x, y, z)

        # Projected target: pyproj returns canonical (east, north) with
        # always_xy=True. The z it returns is unmodified — we ignore it and
        # compute z ourselves from the depth.
        east, north, _z_unused = self._transformer.transform(longitude, latitude, -depth_m)
        # Canonical NED triple (north, east, down):
        ned = (north, east, depth_m)
        return Vec3(*_permute_ned(ned, self.axis))

    def to_mcs(
        self,
        strike: float, dip: float, rake: float,
        u1: float, u2: float, u3: float,
    ) -> Vec3:
        """Rotate a fault-local vector into the mesh CS.

        ``(u1, u2, u3)`` are the components of the vector in the fault-local
        frame: along strike, along dip (down-dip in the fault plane), along
        the fault normal. The rotation goes via canonical NED and is then
        permuted/signed according to ``self.axis``.
        """
        ned = _strike_dip_rake_to_ned(strike, dip, rake, u1, u2, u3)
        return Vec3(*_permute_ned(ned, self.axis))

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------

    def _maybe_warn_units(self) -> None:
        # axis_info has a `unit_name` per axis. We expect 'metre' across the
        # board; anything else means the user mixed units, and downstream
        # area/sliprate conversions (cm^2 → m^2, cm/s → m/s) will produce
        # mismatched output. The C++ tool emits the same warning.
        for axis in self.target_crs.axis_info:
            if axis.unit_name and axis.unit_name not in ("metre", "meter"):
                _log.warning(
                    "Mesh coordinate system does not use metre as length unit "
                    "(axis %s is in %s). NRF area/sliprate are still written in "
                    "m^2 / m/s — expect inconsistencies.",
                    axis.abbrev or axis.direction, axis.unit_name,
                )
                return


# ---------------------------------------------------------------------------
# free functions: rotation and permutation
# ---------------------------------------------------------------------------

def _strike_dip_rake_to_ned(
    strike: float, dip: float, rake: float,
    u1: float, u2: float, u3: float,
) -> tuple[float, float, float]:
    """Rotate ``(u1, u2, u3)`` from the fault-local frame into NED.

    The rotation matrix is the standard one for the SRF convention:

    * ``u1`` is along strike (in the fault plane, in the strike direction)
    * ``u2`` is in the fault plane, perpendicular to strike, pointing down-dip
    * ``u3`` is the fault normal

    The output is in the canonical (north, east, down) frame.

    The expression is verbatim from the C++ ``Map::toMCS`` (before its
    ``adjustAxes`` step), kept in this explicit form so the source-to-source
    correspondence is easy to verify.
    """
    s = math.radians(strike)
    d = math.radians(dip)
    r = math.radians(rake)
    sin_s, cos_s = math.sin(s), math.cos(s)
    sin_d, cos_d = math.sin(d), math.cos(d)
    sin_r, cos_r = math.sin(r), math.cos(r)

    n = (u1 * (sin_r * sin_s * cos_d + cos_r * cos_s)
         + u2 * (-sin_r * cos_s + sin_s * cos_d * cos_r)
         - u3 * sin_d * sin_s)
    e = (u1 * (-sin_r * cos_d * cos_s + sin_s * cos_r)
         + u2 * (-sin_r * sin_s - cos_d * cos_r * cos_s)
         + u3 * sin_d * cos_s)
    d_comp = (-u1 * sin_d * sin_r
              - u2 * sin_d * cos_r
              - u3 * cos_d)
    return n, e, d_comp


def _permute_ned(ned: tuple[float, float, float], axis: str) -> tuple[float, float, float]:
    """Permute and sign a canonical NED triple according to a PROJ ``+axis=`` string.

    Each character of ``axis`` (length 3, characters in ``enwsud``) selects
    which canonical component, with which sign, ends up in the output
    position. For example ``axis='enu'`` produces ``(east, north, -down)``.
    """
    n, e, d = ned
    out = [0.0, 0.0, 0.0]
    for i, ch in enumerate(axis):
        if ch == "n":   out[i] = n
        elif ch == "s": out[i] = -n
        elif ch == "e": out[i] = e
        elif ch == "w": out[i] = -e
        elif ch == "d": out[i] = d
        elif ch == "u": out[i] = -d
        else:
            # Caller validated up-front, but keep this branch for safety.
            raise ValueError(f"Invalid axis character: {ch!r}")
    return out[0], out[1], out[2]


# ---------------------------------------------------------------------------
# proj-string helpers
# ---------------------------------------------------------------------------

_AXIS_RE = re.compile(r"\+axis=(\S+)")


def _parse_axis(proj_string: str) -> str:
    """Extract the ``+axis=`` value from a PROJ string, defaulting to ``"enu"``."""
    m = _AXIS_RE.search(proj_string)
    return m.group(1).lower() if m else _DEFAULT_AXIS
