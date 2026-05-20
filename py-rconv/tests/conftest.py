"""Shared pytest fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

# Make tests/data importable as a module and src/ importable as the rconv pkg
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "data"))
sys.path.insert(0, str(HERE.parent / "src"))

import pytest  # noqa: E402  (must come after sys.path edits)

import synthetic  # noqa: E402  (from tests/data/synthetic.py)


@pytest.fixture
def srf_v10_minimal(tmp_path: Path) -> Path:
    return synthetic.write_srf(synthetic.SRF_V10_MINIMAL, tmp_path / "v10.srf")


@pytest.fixture
def srf_v20_two_dir(tmp_path: Path) -> Path:
    return synthetic.write_srf(synthetic.SRF_V20_TWO_DIR, tmp_path / "v20.srf")


@pytest.fixture
def srf_v10_for_slip_integration(tmp_path: Path) -> Path:
    return synthetic.write_srf(
        synthetic.SRF_V10_FOR_SLIP_INTEGRATION, tmp_path / "slip.srf",
    )


@pytest.fixture
def northridge_srf() -> Path:
    """Real-world Northridge SRF file, if available on this system.

    Skips the test if the file is not present; it lives in the SeisSol
    Examples repository.
    """
    candidates = [
        Path("/home/claude/Examples/Northridge/nr6.70-s0000-h0000.txt"),
        Path(__file__).parent / "data" / "northridge.srf",
    ]
    for p in candidates:
        if p.exists():
            return p
    pytest.skip("Northridge SRF file not available; clone SeisSol/Examples to enable.")
