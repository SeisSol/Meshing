"""Tests for :mod:`rconv.nrf`: round-trip writes, unit conversions, structure."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

from rconv import CoordinateMapping, parse_srf, read_nrf, write_nrf
from rconv.srf import PointSource


@pytest.fixture
def mapping() -> CoordinateMapping:
    return CoordinateMapping(
        "+proj=tmerc +datum=WGS84 +k=0.9996 "
        "+lon_0=-118.5150 +lat_0=34.3440 +units=m +axis=enu"
    )


@pytest.fixture
def single_source() -> PointSource:
    # Hand-crafted source where every quantity is non-zero so unit conversions
    # are easy to check.
    return PointSource(
        longitude=-118.5,
        latitude=34.3,
        depth=5.0,                    # km
        strike=0.0, dip=90.0, rake=0.0,
        area=1.0e8,                   # cm^2 → expect 1e4 m^2
        tinit=1.0, dt=0.5,
        shear_modulus=3.0e10,         # g/(cm*s^2) → expect 3e9 Pa
        slip_rates=([100.0, 50.0], [], []),  # cm/s → expect [1.0, 0.5] m/s
    )


class TestRefuseEmpty:
    def test_empty_sources_raises(self, tmp_path: Path, mapping: CoordinateMapping) -> None:
        with pytest.raises(ValueError, match="empty source list"):
            write_nrf(tmp_path / "out.nrf", [], mapping)


class TestUnitConversion:
    def test_area_cm2_to_m2(
        self, tmp_path: Path, mapping: CoordinateMapping, single_source: PointSource,
    ) -> None:
        out = tmp_path / "out.nrf"
        write_nrf(out, [single_source], mapping)
        data = read_nrf(out)
        # area: 1e8 cm^2 = 1e4 m^2
        assert data["subfaults"][0]["area"] == pytest.approx(1.0e4)

    def test_shear_modulus_cgs_to_pa(
        self, tmp_path: Path, mapping: CoordinateMapping, single_source: PointSource,
    ) -> None:
        out = tmp_path / "out.nrf"
        write_nrf(out, [single_source], mapping)
        data = read_nrf(out)
        # mu: 3e10 g/(cm*s^2) = 3e9 Pa
        assert data["subfaults"][0]["mu"] == pytest.approx(3.0e9)

    def test_sliprate_cm_s_to_m_s(
        self, tmp_path: Path, mapping: CoordinateMapping, single_source: PointSource,
    ) -> None:
        out = tmp_path / "out.nrf"
        write_nrf(out, [single_source], mapping)
        data = read_nrf(out)
        # [100, 50] cm/s → [1.0, 0.5] m/s
        assert list(data["sliprates1"]) == pytest.approx([1.0, 0.5])

    def test_depth_km_to_m(
        self, tmp_path: Path, mapping: CoordinateMapping, single_source: PointSource,
    ) -> None:
        out = tmp_path / "out.nrf"
        write_nrf(out, [single_source], mapping)
        data = read_nrf(out)
        # Depth 5.0 km, axis=enu (z=up) → z = -5000 m
        assert data["centres"][0]["z"] == pytest.approx(-5000.0)


class TestSroffsets:
    """The ``sroffsets`` array must be the cumulative count per direction."""

    def test_cumulative_per_direction(
        self, tmp_path: Path, mapping: CoordinateMapping,
    ) -> None:
        sources = [
            PointSource(0, 0, 1, 0, 90, 0, 1e8, 0.0, 0.1, 0.0, ([1, 2], [3], [])),
            PointSource(0, 0, 1, 0, 90, 0, 1e8, 0.0, 0.1, 0.0, ([4, 5, 6], [], [7, 8])),
            PointSource(0, 0, 1, 0, 90, 0, 1e8, 0.0, 0.1, 0.0, ([9], [], [])),
        ]
        out = tmp_path / "out.nrf"
        write_nrf(out, sources, mapping)
        data = read_nrf(out)
        # Direction 0: 2 + 3 + 1 = 6 samples → offsets 0, 2, 5, 6
        # Direction 1:     1 + 0 + 0 = 1 samples → offsets 0, 1, 1, 1
        # Direction 2:     0 + 2 + 0 = 2 samples → offsets 0, 0, 2, 2
        expected = np.array([
            [0, 0, 0],
            [2, 1, 0],
            [5, 1, 2],
            [6, 1, 2],
        ], dtype=np.uint32)
        assert np.array_equal(data["sroffsets"], expected)

    def test_final_offset_equals_sliprate_length(
        self, tmp_path: Path, mapping: CoordinateMapping, single_source: PointSource,
    ) -> None:
        out = tmp_path / "out.nrf"
        write_nrf(out, [single_source], mapping)
        data = read_nrf(out)
        # For n sources, sroffsets[n, d] must equal sliprates{d+1}.size.
        n = data["centres"].size
        assert data["sroffsets"][n, 0] == data["sliprates1"].size
        assert data["sroffsets"][n, 1] == data["sliprates2"].size
        assert data["sroffsets"][n, 2] == data["sliprates3"].size


class TestNormalizeOnset:
    def test_subtracts_minimum_tinit(
        self, tmp_path: Path, mapping: CoordinateMapping,
    ) -> None:
        # PointSource positional order:
        # lon, lat, depth, strike, dip, rake, area, tinit, dt, shear_modulus, slip_rates
        sources = [
            PointSource(0, 0, 1, 0, 90, 0, 1e8, 5.0, 0.1, 0.0, ([1.0], [], [])),
            PointSource(0, 0, 1, 0, 90, 0, 1e8, 7.5, 0.1, 0.0, ([1.0], [], [])),
            PointSource(0, 0, 1, 0, 90, 0, 1e8, 3.0, 0.1, 0.0, ([1.0], [], [])),  # min
        ]
        out = tmp_path / "out.nrf"
        write_nrf(out, sources, mapping, normalize_onset=True)
        data = read_nrf(out)
        assert data["subfaults"][0]["tinit"] == pytest.approx(2.0)
        assert data["subfaults"][1]["tinit"] == pytest.approx(4.5)
        assert data["subfaults"][2]["tinit"] == pytest.approx(0.0)

    def test_off_by_default(
        self, tmp_path: Path, mapping: CoordinateMapping, single_source: PointSource,
    ) -> None:
        out = tmp_path / "out.nrf"
        write_nrf(out, [single_source], mapping)
        data = read_nrf(out)
        # Default — onsets pass through unmodified.
        assert data["subfaults"][0]["tinit"] == pytest.approx(1.0)


class TestStructure:
    """Smoke tests for the on-disk file structure via ncdump."""

    def test_compound_types_present(
        self, tmp_path: Path, mapping: CoordinateMapping, single_source: PointSource,
    ) -> None:
        if not _have_ncdump():
            pytest.skip("ncdump not on PATH")
        out = tmp_path / "out.nrf"
        write_nrf(out, [single_source], mapping)
        header = subprocess.check_output(
            ["ncdump", "-h", str(out)], text=True,
        )
        assert "compound Vector3" in header
        assert "compound Subfault" in header
        # Expected vars
        for var in ("centres", "subfaults", "sroffsets",
                    "sliprates1", "sliprates2", "sliprates3"):
            assert var in header


class TestRealistic:
    def test_northridge_round_trip(
        self, tmp_path: Path, northridge_srf: Path, mapping: CoordinateMapping,
    ) -> None:
        sources = parse_srf(northridge_srf)
        out = tmp_path / "northridge.nrf"
        write_nrf(out, sources, mapping)
        data = read_nrf(out)
        assert data["centres"].size == 474
        assert data["sroffsets"].shape == (475, 3)
        # All depths positive in SRF → all z negative in ENU output.
        assert (data["centres"]["z"] <= 0).all()
        # Every source's tan1 / tan2 / normal should be a unit vector triple.
        for sf in data["subfaults"]:
            for name in ("tan1", "tan2", "normal"):
                v = sf[name]
                norm = (v["x"] ** 2 + v["y"] ** 2 + v["z"] ** 2) ** 0.5
                assert norm == pytest.approx(1.0, abs=1e-10)


def _have_ncdump() -> bool:
    try:
        # ncdump prints usage and exits 0 when given no args.
        subprocess.run(
            ["ncdump"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            check=False,
        )
        return True
    except FileNotFoundError:
        return False
