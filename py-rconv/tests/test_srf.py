"""Tests for the SRF parser."""

from __future__ import annotations

from pathlib import Path

import pytest

from rconv.srf import PointSource, SRFParseError, parse_srf


class TestVersionDispatch:
    def test_unsupported_version_raises(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.srf"
        bad.write_text("0.5\nPOINTS 0\n")
        with pytest.raises(SRFParseError, match="unsupported SRF version"):
            parse_srf(bad)

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty.srf"
        empty.write_text("")
        with pytest.raises(SRFParseError, match="empty file"):
            parse_srf(empty)


class TestV10:
    def test_minimal_parse_drops_empty_sources(self, srf_v10_minimal: Path) -> None:
        # The minimal fixture has 2 sources, the second of which carries no
        # slip-rate samples in any direction → must be dropped.
        sources = parse_srf(srf_v10_minimal)
        assert len(sources) == 1
        s = sources[0]
        assert s.longitude == -118.5
        assert s.latitude == 34.3
        assert s.depth == 5.0

    def test_v10_has_no_shear_modulus(self, srf_v10_minimal: Path) -> None:
        # v1.0 doesn't carry vs/den, so shear modulus stays 0 (unknown).
        s = parse_srf(srf_v10_minimal)[0]
        assert s.shear_modulus == 0.0

    def test_v10_samples_have_expected_counts(self, srf_v10_minimal: Path) -> None:
        s = parse_srf(srf_v10_minimal)[0]
        assert len(s.slip_rates[0]) == 4
        assert len(s.slip_rates[1]) == 0
        assert len(s.slip_rates[2]) == 0

    def test_v10_sample_values(self, srf_v10_minimal: Path) -> None:
        s = parse_srf(srf_v10_minimal)[0]
        assert s.slip_rates[0] == [0.0, 1.0, 2.0, 0.0]


class TestV20:
    def test_v20_shear_modulus_from_vs_den(self, srf_v20_two_dir: Path) -> None:
        # vs=350000 cm/s, den=2.7 g/cm^3 → mu = vs^2 * den
        s = parse_srf(srf_v20_two_dir)[0]
        expected = 350000.0 ** 2 * 2.7
        assert s.shear_modulus == pytest.approx(expected)

    def test_v20_plane_block_skipped(self, srf_v20_two_dir: Path) -> None:
        # If the PLANE block were not consumed correctly the subsequent POINTS
        # parse would fail or produce garbage. Reaching here with the right
        # number of sources is the assertion.
        sources = parse_srf(srf_v20_two_dir)
        assert len(sources) == 1
        # Comment line was also handled correctly.
        assert sources[0].latitude == 34.3

    def test_v20_two_direction_samples(self, srf_v20_two_dir: Path) -> None:
        s = parse_srf(srf_v20_two_dir)[0]
        assert s.slip_rates[0] == [1.0, 2.0, 1.5]
        assert s.slip_rates[1] == [0.5, 0.25]
        assert s.slip_rates[2] == []


class TestV20EdgeCases:
    def test_v20_zero_vs_or_den_gives_unknown_mu(self, tmp_path: Path) -> None:
        # The C++ tool treats vs=0 or den=0 as "unknown" → mu = 0.
        srf = tmp_path / "v20_zero.srf"
        srf.write_text(
            "2.0\n"
            "POINTS 1\n"
            "0.0 0.0 1.0 0.0 90.0 1.0e4 0.0 0.01 0.0 2.7\n"  # vs=0
            "0.0 1.0 1 0.0 0 0.0 0\n"
            "1.0\n"
        )
        s = parse_srf(srf)[0]
        assert s.shear_modulus == 0.0


class TestRealistic:
    def test_northridge_loads(self, northridge_srf: Path) -> None:
        # The Northridge SRF has 500 raw point sources, 26 of which are
        # empty placeholders the parser must drop.
        sources = parse_srf(northridge_srf)
        assert len(sources) == 474

    def test_northridge_first_source(self, northridge_srf: Path) -> None:
        s = parse_srf(northridge_srf)[0]
        # Values straight from the SRF; if any of these change the parser is
        # mis-reading the file.
        assert s.longitude == pytest.approx(-118.6049)
        assert s.latitude == pytest.approx(34.3864)
        assert s.depth == pytest.approx(5.3214)
        assert s.strike == pytest.approx(122.0)
        assert s.dip == pytest.approx(40.0)
        assert s.rake == pytest.approx(90.0)
        assert len(s.slip_rates[0]) == 16
