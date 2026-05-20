"""Tests for :mod:`rconv.xmf`."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from rconv import CoordinateMapping, parse_srf, write_xmf
from rconv.xmf import _slip_path_length_m


@pytest.fixture
def mapping() -> CoordinateMapping:
    return CoordinateMapping("+proj=geocent +datum=WGS84 +units=m +no_defs")


class TestStructure:
    def test_is_valid_xml(
        self, tmp_path: Path, srf_v10_for_slip_integration: Path, mapping: CoordinateMapping,
    ) -> None:
        sources = parse_srf(srf_v10_for_slip_integration)
        out = tmp_path / "out.xmf"
        write_xmf(out, sources, mapping)
        # ET parses xmf even with the DOCTYPE present.
        tree = ET.parse(out)
        root = tree.getroot()
        assert root.tag == "Xdmf"

    def test_polyvertex_count_matches(
        self, tmp_path: Path, northridge_srf: Path, mapping: CoordinateMapping,
    ) -> None:
        sources = parse_srf(northridge_srf)
        out = tmp_path / "nr.xmf"
        write_xmf(out, sources, mapping)
        tree = ET.parse(out)
        topo = tree.getroot().find(".//Topology")
        assert topo is not None
        assert topo.get("Type") == "Polyvertex"
        assert int(topo.get("NumberOfElements", "0")) == len(sources)

    def test_has_slip_attribute(
        self, tmp_path: Path, srf_v10_for_slip_integration: Path, mapping: CoordinateMapping,
    ) -> None:
        sources = parse_srf(srf_v10_for_slip_integration)
        out = tmp_path / "out.xmf"
        write_xmf(out, sources, mapping)
        tree = ET.parse(out)
        attr = tree.getroot().find(".//Attribute")
        assert attr is not None
        assert attr.get("Name") == "Slip path length (m)"
        assert attr.get("Center") == "Node"


class TestSlipPathLength:
    """The XDMF 'Slip path length' attribute is the trapezoidal integral of
    the slip-rate magnitude over time, then converted from cm to m."""

    def test_trapezoid_known_value(
        self, srf_v10_for_slip_integration: Path,
    ) -> None:
        sources = parse_srf(srf_v10_for_slip_integration)
        s = sources[0]
        # The fixture has u1 = [1, 2, 3, 2, 1] cm/s, dt = 0.5 s.
        # Trapezoidal integral with dx=0.5 on [1,2,3,2,1] = 4.0 cm = 0.04 m.
        expected_cm = float(np.trapezoid([1.0, 2.0, 3.0, 2.0, 1.0], dx=0.5))
        assert expected_cm == pytest.approx(4.0)
        assert _slip_path_length_m(s) == pytest.approx(0.04)

    def test_empty_source_yields_zero(self) -> None:
        from rconv.srf import PointSource
        empty = PointSource(
            longitude=0, latitude=0, depth=0,
            strike=0, dip=90, rake=0,
            area=1e8, tinit=0, dt=0.1, shear_modulus=0,
            slip_rates=([], [], []),
        )
        assert _slip_path_length_m(empty) == 0.0

    def test_multi_direction_magnitude(self) -> None:
        # If u1=[3,3] and u3=[4,4] cm/s, the vector magnitude is 5 at every
        # sample. Trapezoid with dx=1.0 over 2 samples = 5*1.0 = 5 cm = 0.05 m.
        from rconv.srf import PointSource
        s = PointSource(
            longitude=0, latitude=0, depth=0,
            strike=0, dip=90, rake=0,
            area=1e8, tinit=0, dt=1.0, shear_modulus=0,
            slip_rates=([3.0, 3.0], [], [4.0, 4.0]),
        )
        assert _slip_path_length_m(s) == pytest.approx(0.05)
