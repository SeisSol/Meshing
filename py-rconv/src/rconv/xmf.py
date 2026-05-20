"""Writer for XDMF/XMF visualisation output.

Produces an XML file that ParaView (and any other XDMF-aware viewer) can open
to render the point cloud of SRF sources, decorated with the integrated slip
path length per source.

The output mirrors the C++ ``rconv``'s ``writeXMF`` exactly in structure (a
``Polyvertex`` topology with inline XYZ geometry and a single ``Slip path
length (m)`` node attribute), modulo whitespace.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

from .mapping import CoordinateMapping
from .srf import PointSource

# cm/s → m/s on slip path: integration is in (cm/s)·s = cm, then /100 → m.
_CM_TO_M = 1.0e-2


def write_xmf(
    path: str | Path,
    sources: Sequence[PointSource],
    mapping: CoordinateMapping,
) -> None:
    """Write an XDMF/XMF visualisation file for a list of point sources.

    Parameters
    ----------
    path:
        Destination file path. Overwritten if it exists.
    sources:
        Point sources, as produced by :func:`rconv.srf.parse_srf`.
    mapping:
        Coordinate mapping used to project source positions for display.
        Typically a geocentric CRS for global views or the same MCS as for the
        NRF for mesh-aligned visualisation.
    """
    path = Path(path)
    n = len(sources)

    with path.open("w") as out:
        out.write('<?xml version="1.0" ?>\n')
        out.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        out.write("<Xdmf>\n")
        out.write("  <Domain>\n")
        out.write(f'    <Topology Type="Polyvertex" NumberOfElements="{n}"/>\n')
        out.write('    <Geometry Type="XYZ">\n')
        out.write(
            f'      <DataItem Dimensions="{n} 3" NumberType="Float" Precision="8" Format="XML">\n'
        )
        for src in sources:
            x, y, z = mapping.to_mesh(src.longitude, src.latitude, src.depth).as_tuple()
            out.write(f"{x} {y} {z}\n")
        out.write("      </DataItem>\n")
        out.write("    </Geometry>\n")
        out.write('    <Grid Name="Fault" GridType="Uniform">\n')
        out.write('      <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n')
        out.write('      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n')
        _write_slip_attribute(out, sources)
        out.write("    </Grid>\n")
        out.write("  </Domain>\n")
        out.write("</Xdmf>\n")


def _write_slip_attribute(out, sources: Sequence[PointSource]) -> None:
    """Write the per-source 'Slip path length (m)' attribute.

    Path length is the time-integral of the magnitude of the (vector) slip
    rate, evaluated with the trapezoidal rule on the SRF sample grid. The C++
    ``rconv`` does the same thing.
    """
    n = len(sources)
    out.write('      <Attribute Name="Slip path length (m)" Center="Node">\n')
    out.write(
        f'        <DataItem Dimensions="{n}" DataType="Float" Precision="8" Format="XML">\n'
    )
    for src in sources:
        out.write(f"{_slip_path_length_m(src)}\n")
    out.write("        </DataItem>\n")
    out.write("      </Attribute>\n")


def _slip_path_length_m(src: PointSource) -> float:
    """Trapezoidal integral of |slip rate vector| over time, in metres.

    Samples are padded with zeros if the per-direction sample counts differ
    (which is permitted by SRF but uncommon in practice).
    """
    max_steps = max((len(r) for r in src.slip_rates), default=0)
    if max_steps == 0:
        return 0.0
    # Pack into an (max_steps, 3) array, zero-padding short components.
    sr = np.zeros((max_steps, 3), dtype=np.float64)
    for d, samples in enumerate(src.slip_rates):
        if samples:
            sr[: len(samples), d] = samples
    magnitudes = np.linalg.norm(sr, axis=1)  # cm/s, per sample
    # Trapezoidal rule with uniform spacing dt; numpy's trapezoid is exactly
    # what the C++ implementation does pointwise.
    slip_cm = float(np.trapezoid(magnitudes, dx=src.dt))
    return slip_cm * _CM_TO_M
