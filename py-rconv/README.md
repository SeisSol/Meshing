# rconv-py

A pure-Python rewrite of SeisSol's `rconv` tool.

`rconv` converts kinematic rupture models from the [Standard Rupture
Format](https://strike.scec.org/scecpedia/Standard_Rupture_Format) (SRF) to
SeisSol's NetCDF Rupture Format (NRF). It can also emit an XDMF file for
visualisation in ParaView.

The original `rconv` is a C++ tool that depends on a deprecated version of
PROJ.4 and on the SeisSol build system. This package replaces it with a small,
self-contained Python module that uses modern `pyproj` and `netCDF4` — nothing
else.

## Installation

```sh
pip install .
```

Dependencies (installed automatically):

* `numpy >= 1.24`
* `pyproj >= 3.4`
* `netCDF4 >= 1.6`

## Quick start

```sh
# Convert SRF → NRF for SeisSol input.
rconv to-nrf input.srf output.nrf \
    --mcs "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=-118.5150 +lat_0=34.3440 +axis=enu"

# Produce an XDMF visualisation file for ParaView.
rconv to-xdmf input.srf visualization.xmf
```

Equivalent invocation as a module (no entry-point installation needed):

```sh
python -m rconv to-nrf input.srf output.nrf --mcs "..."
```

## CLI

Two subcommands cover the two distinct conversions:

```
rconv to-nrf INPUT OUTPUT --mcs PROJ_STRING [--normalize-onset] [--source-proj PROJ_STRING]
rconv to-xdmf INPUT OUTPUT [--vcs PROJ_STRING] [--source-proj PROJ_STRING]
```

* `--mcs` — PROJ string for the *mesh* coordinate system. **Required** for NRF
  output. Determines both how lon/lat/depth becomes (x, y, z) and how
  slip-vector components are oriented (via the `+axis=` option).
* `--vcs` — PROJ string for the *visualisation* coordinate system. Defaults to
  `+proj=geocent +datum=WGS84 +units=m +no_defs` (a geocentric globe view).
* `--normalize-onset` — subtract the minimum onset time across all sources.
  Useful when the SRF's absolute clock is meaningless.
* `--source-proj` — override the source CRS (default: WGS84 longlat). Useful
  only when the SRF has already-projected coordinates.

## Programmatic API

```python
from rconv import parse_srf, CoordinateMapping, write_nrf, write_xmf

sources = parse_srf("northridge.srf")
mapping = CoordinateMapping(
    "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=-118.515 +lat_0=34.344 +axis=ned"
)
write_nrf("northridge.nrf", sources, mapping)
write_xmf("northridge.xmf", sources, mapping)
```

## What `+axis=` means

The `+axis=` value in your PROJ string controls how a canonical
(north, east, down) triple — the natural frame of SRF coordinates — is
permuted and signed into the mesh CS. It is three characters from `enwsud`:

* `+axis=enu` — x = east, y = north, z = up (PROJ's default if you omit
  `+axis=` entirely).
* `+axis=ned` — x = north, y = east, z = down.
* `+axis=seu` — x = south, y = east, z = up.

The same string also controls how slip vectors (along strike, in fault plane,
along normal) are oriented in your mesh. If your moment tensors look wrong,
double-check `+axis=` first.

## Differences from the C++ tool

* **CLI**: two focused subcommands instead of five short flags.
* **PROJ**: uses modern `pyproj` (which wraps PROJ 8/9), not the legacy
  PROJ.4 API the C++ tool depended on.
* **Output**: NRF output is functionally identical and consumed by SeisSol's
  NRF reader without modification. Two cosmetic differences:
  - The `Subfault_units` compound type is replaced by a plain string
    attribute on `subfaults`. SeisSol's reader never inspects either form,
    so this is invisible to downstream use; only `ncdump` output differs.
  - Slip-rate dimensions of size 0 are written as UNLIMITED rather than
    fixed-size-0. Again, the reader doesn't care.
* **Examples**: ships with a test suite that round-trips the Northridge SRF
  from `SeisSol/Examples`.

## Development

Run the tests:

```sh
pip install -e .[dev]
pytest
```

The Northridge integration tests are skipped automatically if the SRF file
isn't available at `~/Examples/Northridge/`.
