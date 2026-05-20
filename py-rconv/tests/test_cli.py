"""Tests for the command-line interface."""

from __future__ import annotations

from pathlib import Path

import pytest

from rconv.cli import main
from rconv.nrf import read_nrf


class TestArgparse:
    def test_no_command_exits_nonzero(self, capsys: pytest.CaptureFixture[str]) -> None:
        # argparse "required=True" on the subparsers slot makes invocation
        # without a subcommand fail with exit code 2.
        with pytest.raises(SystemExit) as exc:
            main([])
        assert exc.value.code == 2

    def test_unknown_subcommand_fails(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["frobnicate", "a", "b"])
        assert exc.value.code == 2

    def test_to_nrf_requires_mcs(self, tmp_path: Path) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["to-nrf", str(tmp_path / "x.srf"), str(tmp_path / "y.nrf")])
        assert exc.value.code == 2


class TestEndToEnd:
    def test_to_nrf_writes_file(
        self, tmp_path: Path, srf_v10_minimal: Path,
    ) -> None:
        out = tmp_path / "out.nrf"
        code = main([
            "to-nrf", str(srf_v10_minimal), str(out),
            "--mcs", "+proj=tmerc +datum=WGS84 +k=0.9996 "
                     "+lon_0=-118.5150 +lat_0=34.3440 +units=m +axis=enu",
        ])
        assert code == 0
        assert out.exists()
        data = read_nrf(out)
        # The fixture's first source is kept, the empty one is dropped.
        assert data["centres"].size == 1

    def test_to_nrf_normalize_onset(
        self, tmp_path: Path, srf_v10_minimal: Path,
    ) -> None:
        out = tmp_path / "out.nrf"
        code = main([
            "to-nrf", str(srf_v10_minimal), str(out),
            "--mcs", "+proj=tmerc +datum=WGS84 +k=0.9996 "
                     "+lon_0=-118.5150 +lat_0=34.3440 +units=m +axis=enu",
            "--normalize-onset",
        ])
        assert code == 0
        data = read_nrf(out)
        # Single source survives → min tinit = its own tinit → result 0.
        assert data["subfaults"][0]["tinit"] == pytest.approx(0.0)

    def test_to_xdmf_writes_file(
        self, tmp_path: Path, srf_v10_minimal: Path,
    ) -> None:
        out = tmp_path / "out.xmf"
        code = main(["to-xdmf", str(srf_v10_minimal), str(out)])
        assert code == 0
        assert out.exists()
        text = out.read_text()
        assert "<Xdmf>" in text
        assert "Polyvertex" in text

    def test_to_xdmf_default_vcs_is_geocentric(
        self, tmp_path: Path, srf_v10_minimal: Path,
    ) -> None:
        # The default VCS is geocentric; result coordinates should have
        # magnitude ~ Earth radius.
        import re
        out = tmp_path / "out.xmf"
        code = main(["to-xdmf", str(srf_v10_minimal), str(out)])
        assert code == 0
        # First float triple after the geometry DataItem opening.
        m = re.search(
            r'Geometry.*?XML">\s*([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)',
            out.read_text(), re.DOTALL,
        )
        assert m is not None
        x, y, z = map(float, m.groups())
        r = (x * x + y * y + z * z) ** 0.5
        assert 6.3e6 < r < 6.4e6


class TestEmptyInput:
    def test_no_sources_returns_error(
        self, tmp_path: Path,
    ) -> None:
        # An SRF whose sole point source has no samples in any direction
        # parses to an empty list and should be a CLI error rather than an
        # unhelpful crash.
        srf = tmp_path / "empty.srf"
        srf.write_text(
            "1.0\n"
            "POINTS 1\n"
            "0 0 1 0 90 1e8 0 0.1\n"
            "0 0 0 0 0 0 0\n"
        )
        code = main([
            "to-nrf", str(srf), str(tmp_path / "out.nrf"),
            "--mcs", "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=0 +lat_0=0 +units=m",
        ])
        assert code == 1
