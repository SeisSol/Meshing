"""Parser for the Standard Rupture Format (SRF).

Reference: https://strike.scec.org/scecpedia/Standard_Rupture_Format

This module reads SRF v1.0 and v2.0 files. Both versions describe a kinematic
rupture model as a collection of point sources, each with a per-source time
history of slip rate in three orthogonal directions of the fault-local frame.

The parser is tokenizer-based (rather than line-based), mirroring the C++
``rconv`` behaviour: SRF files mix whitespace and newlines as separators within
both the per-point header lines and the slip-rate sample arrays.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from io import TextIOBase
from pathlib import Path
from typing import IO, Iterator

_log = logging.getLogger(__name__)


@dataclass
class PointSource:
    """A single SRF point source.

    Field units follow the SRF specification verbatim — converting them to SI
    is the responsibility of downstream code (e.g. :func:`rconv.nrf.write_nrf`).
    """

    longitude: float        #: degrees
    latitude: float         #: degrees
    depth: float            #: km, positive downwards
    strike: float           #: degrees
    dip: float              #: degrees
    rake: float             #: degrees
    area: float             #: cm^2
    tinit: float            #: s, onset time of the source-time function
    dt: float               #: s, sample spacing of the slip-rate function
    shear_modulus: float    #: g/(cm*s^2); 0.0 if unknown (v1.0 always; v2.0 if vs or den missing)
    slip_rates: tuple[list[float], list[float], list[float]] = field(
        default_factory=lambda: ([], [], [])
    )                       #: cm/s, samples per fault-local direction (u1, u2, u3)


class SRFParseError(ValueError):
    """Raised when an SRF file cannot be parsed."""


def parse_srf(path: str | Path) -> list[PointSource]:
    """Parse an SRF file and return the list of (non-empty) point sources.

    Point sources whose slip-rate samples are empty in all three directions are
    dropped — they carry no information and would otherwise cause divide-by-zero
    issues downstream. The number removed is logged at INFO level.

    Parameters
    ----------
    path:
        Path to the SRF file.

    Returns
    -------
    list[PointSource]
        The parsed point sources, in file order.

    Raises
    ------
    SRFParseError
        If the file has an unsupported version, an unknown block type, or is
        malformed.
    """
    path = Path(path)
    with path.open("r") as fh:
        return _parse_stream(fh, source=str(path))


def _parse_stream(stream: IO[str] | TextIOBase, *, source: str = "<stream>") -> list[PointSource]:
    tokens = _tokenize(stream)

    # First token: version. Accept floats within tolerance, as the C++ tool
    # parses them with `in >> version` and then compares to 1.0 / 2.0 exactly.
    try:
        version = float(next(tokens))
    except StopIteration as exc:
        raise SRFParseError(f"{source}: empty file") from exc

    if abs(version - 1.0) < 1e-6:
        parse_point = _parse_point_v10
    elif abs(version - 2.0) < 1e-6:
        parse_point = _parse_point_v20
    else:
        raise SRFParseError(f"{source}: unsupported SRF version {version}")

    sources: list[PointSource] = []

    for block in tokens:
        if block == "POINTS":
            n_points = int(next(tokens))
            for _ in range(n_points):
                sources.append(parse_point(tokens))
        elif block == "PLANE":
            # The PLANE block has 2 header lines per plane in v1.0/v2.0. The
            # C++ tool skips them entirely; we do the same since the per-point
            # records following carry all the information we actually need.
            n_planes = int(next(tokens))
            _skip_plane(tokens, n_planes)
        elif block.startswith("#"):
            # Comments span until end of line; the tokenizer already strips them
            # (see `_tokenize`), so this branch shouldn't normally fire — kept
            # as a defensive fallback in case a `#` appears mid-stream.
            continue
        else:
            raise SRFParseError(f"{source}: unknown block keyword {block!r}")

    # Drop point sources that contain no samples whatsoever.
    keep = [s for s in sources if any(s.slip_rates)]
    removed = len(sources) - len(keep)
    if removed:
        _log.info("Removed %d point source(s) with no samples in any direction.", removed)
    return keep


# ---------------------------------------------------------------------------
# tokenizer
# ---------------------------------------------------------------------------

def _tokenize(stream: IO[str] | TextIOBase) -> Iterator[str]:
    """Yield whitespace-separated tokens from `stream`, stripping `#` comments.

    The SRF format treats newlines and spaces interchangeably within records,
    so a flat token stream is the natural representation. Comments start with
    `#` and extend to the end of the line, exactly as in the C++ rconv tool.
    """
    for line in stream:
        # Strip everything from the first `#` onward (comment to end-of-line).
        hash_idx = line.find("#")
        if hash_idx != -1:
            line = line[:hash_idx]
        yield from line.split()


# ---------------------------------------------------------------------------
# block readers
# ---------------------------------------------------------------------------

def _skip_plane(tokens: Iterator[str], n_planes: int) -> None:
    """Consume the two header records that follow a `PLANE n` block.

    Each PLANE block contains, per plane, two records of fixed length:

      line 1:  ELON ELAT NSTK NDIP LEN WID
      line 2:  STK DIP DTOP SHYP DHYP

    That's 6 + 5 = 11 tokens per plane. We don't use these for point-source
    conversion (the per-point records following carry their own geometry).
    """
    _consume(tokens, 11 * n_planes)


def _parse_point_v10(tokens: Iterator[str]) -> PointSource:
    # Line 1: LON LAT DEP STK DIP AREA TINIT DT
    lon, lat, dep, stk, dip, area, tinit, dt = _consume_floats(tokens, 8)
    # Line 2: RAKE SLIP1 NT1 SLIP2 NT2 SLIP3 NT3
    rake, _slip1, nt1, _slip2, nt2, _slip3, nt3 = _consume_floats(tokens, 7)
    # v1.0 does not provide rho/vs, so shear modulus is unknown.
    return _read_point(
        tokens,
        lon=lon, lat=lat, dep=dep, stk=stk, dip=dip, rake=rake,
        area=area, tinit=tinit, dt=dt, shear_modulus=0.0,
        nt=(int(nt1), int(nt2), int(nt3)),
    )


def _parse_point_v20(tokens: Iterator[str]) -> PointSource:
    # Line 1: LON LAT DEP STK DIP AREA TINIT DT VS DEN
    lon, lat, dep, stk, dip, area, tinit, dt, vs, den = _consume_floats(tokens, 10)
    # Line 2: RAKE SLIP1 NT1 SLIP2 NT2 SLIP3 NT3
    rake, _slip1, nt1, _slip2, nt2, _slip3, nt3 = _consume_floats(tokens, 7)
    # Shear modulus mu = rho * vs^2. Units: vs in cm/s (per SRF), den in g/cm^3
    # → mu in g/(cm*s^2). The C++ tool uses the same product and silently
    # treats non-positive vs/den as "unknown" (mu = 0).
    shear_modulus = vs * vs * den if vs > 0.0 and den > 0.0 else 0.0
    return _read_point(
        tokens,
        lon=lon, lat=lat, dep=dep, stk=stk, dip=dip, rake=rake,
        area=area, tinit=tinit, dt=dt, shear_modulus=shear_modulus,
        nt=(int(nt1), int(nt2), int(nt3)),
    )


def _read_point(
    tokens: Iterator[str],
    *,
    lon: float, lat: float, dep: float,
    stk: float, dip: float, rake: float,
    area: float, tinit: float, dt: float,
    shear_modulus: float,
    nt: tuple[int, int, int],
) -> PointSource:
    slip_rates = tuple(_consume_floats(tokens, n) for n in nt)  # type: ignore[assignment]
    return PointSource(
        longitude=lon, latitude=lat, depth=dep,
        strike=stk, dip=dip, rake=rake,
        area=area, tinit=tinit, dt=dt,
        shear_modulus=shear_modulus,
        slip_rates=slip_rates,
    )


# ---------------------------------------------------------------------------
# small token helpers
# ---------------------------------------------------------------------------

def _consume(tokens: Iterator[str], n: int) -> list[str]:
    out: list[str] = []
    for _ in range(n):
        try:
            out.append(next(tokens))
        except StopIteration as exc:
            raise SRFParseError(
                f"unexpected end of file (wanted {n} tokens, got {len(out)})"
            ) from exc
    return out


def _consume_floats(tokens: Iterator[str], n: int) -> list[float]:
    raw = _consume(tokens, n)
    try:
        return [float(x) for x in raw]
    except ValueError as exc:
        raise SRFParseError(f"expected {n} floats, got: {raw!r}") from exc
