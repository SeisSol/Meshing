"""Modern command-line interface for ``rconv``.

The classic C++ tool accepts five short flags that mix concerns (some required
together, others mutually independent). Here we split into two focused
sub-commands so the help is self-documenting and shell-completion is happy:

* ``rconv to-nrf INPUT OUTPUT --mcs ...`` — produce a SeisSol input file.
* ``rconv to-xdmf INPUT OUTPUT [--vcs ...]`` — produce a visualisation file.

Both commands accept ``--source-proj`` to override the default input CRS
(WGS84 longlat) on the rare occasion an SRF arrives in something else.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .mapping import CoordinateMapping
from .nrf import write_nrf
from .srf import parse_srf
from .xmf import write_xmf

_DEFAULT_VCS = "+proj=geocent +datum=WGS84 +units=m +no_defs"


def main(argv: list[str] | None = None) -> int:
    """Entry point. Returns a process exit code."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if args.command == "to-nrf":
        return _cmd_to_nrf(args)
    if args.command == "to-xdmf":
        return _cmd_to_xdmf(args)
    parser.print_help()
    return 2


# ---------------------------------------------------------------------------
# argparse plumbing
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rconv",
        description=(
            "Convert SRF (Standard Rupture Format) files to NRF (NetCDF "
            "Rupture Format) for SeisSol, or to XDMF for visualisation."
        ),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging.",
    )
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    nrf = sub.add_parser(
        "to-nrf",
        help="Convert SRF → NRF.",
        description="Convert an SRF file to a NetCDF Rupture Format file for SeisSol.",
    )
    nrf.add_argument("input", type=Path, help="Input SRF file (.srf).")
    nrf.add_argument("output", type=Path, help="Output NRF file (.nrf).")
    nrf.add_argument(
        "--mcs", required=True, metavar="PROJ_STRING",
        help=(
            "PROJ string describing the mesh coordinate system, e.g. "
            '"+proj=utm +zone=10 +datum=WGS84 +units=m +axis=ned +no_defs". '
            "The +axis= option determines how slip vectors are oriented in the "
            "MCS (default: enu)."
        ),
    )
    nrf.add_argument(
        "--normalize-onset", action="store_true",
        help="Subtract the minimum onset time across sources from every source.",
    )
    nrf.add_argument(
        "--source-proj", default=None, metavar="PROJ_STRING",
        help=(
            "Override the source CRS (defaults to WGS84 longlat). Useful only "
            "if your SRF carries already-projected coordinates."
        ),
    )

    xmf = sub.add_parser(
        "to-xdmf",
        help="Convert SRF → XDMF (for ParaView).",
        description="Convert an SRF file to an XDMF visualisation file.",
    )
    xmf.add_argument("input", type=Path, help="Input SRF file (.srf).")
    xmf.add_argument("output", type=Path, help="Output XDMF file (.xmf/.xdmf).")
    xmf.add_argument(
        "--vcs", default=_DEFAULT_VCS, metavar="PROJ_STRING",
        help=(
            "PROJ string for the visualisation coordinate system. "
            f"Defaults to {_DEFAULT_VCS!r} (geocentric — gives a globe view)."
        ),
    )
    xmf.add_argument(
        "--source-proj", default=None, metavar="PROJ_STRING",
        help="Override the source CRS (defaults to WGS84 longlat).",
    )

    return parser


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def _cmd_to_nrf(args: argparse.Namespace) -> int:
    sources = _read_srf(args.input)
    if not sources:
        print(f"error: {args.input} produced no usable point sources", file=sys.stderr)
        return 1
    mapping = CoordinateMapping(args.mcs, source_proj=args.source_proj)
    write_nrf(args.output, sources, mapping, normalize_onset=args.normalize_onset)
    print(f"Wrote {args.output} ({len(sources)} sources).")
    return 0


def _cmd_to_xdmf(args: argparse.Namespace) -> int:
    sources = _read_srf(args.input)
    if not sources:
        print(f"error: {args.input} produced no usable point sources", file=sys.stderr)
        return 1
    mapping = CoordinateMapping(args.vcs, source_proj=args.source_proj)
    write_xmf(args.output, sources, mapping)
    print(f"Wrote {args.output} ({len(sources)} sources).")
    return 0


def _read_srf(path: Path) -> list:
    print(f"Reading {path} ...", end=" ", flush=True)
    sources = parse_srf(path)
    print(f"{len(sources)} sources.")
    return sources


if __name__ == "__main__":
    sys.exit(main())
