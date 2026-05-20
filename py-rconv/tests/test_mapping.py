"""Tests for :mod:`rconv.mapping`."""

from __future__ import annotations

import logging
import math

import pytest

from rconv.mapping import (
    CoordinateMapping,
    Vec3,
    _parse_axis,
    _permute_ned,
    _strike_dip_rake_to_ned,
)


class TestAxisParsing:
    @pytest.mark.parametrize(
        "proj, expected",
        [
            ("+proj=utm +zone=10 +axis=ned", "ned"),
            ("+proj=utm +zone=10", "enu"),                  # default
            ("+axis=SeU +proj=tmerc", "seu"),               # case-insensitive
            ("+proj=geocent +axis=enu +datum=WGS84", "enu"),
        ],
    )
    def test_parse_axis(self, proj: str, expected: str) -> None:
        assert _parse_axis(proj) == expected

    def test_invalid_axis_chars_rejected_by_proj(self) -> None:
        # PROJ itself rejects unrecognised axis characters before our own
        # check fires, but the effect is the same — construction fails.
        with pytest.raises(Exception, match="(?i)axis|proj"):
            CoordinateMapping("+proj=utm +zone=10 +datum=WGS84 +units=m +axis=xyz")

    def test_permute_ned_rejects_bad_char(self) -> None:
        # _permute_ned has its own guard for defence in depth.
        with pytest.raises(ValueError, match="Invalid axis character"):
            _permute_ned((1.0, 2.0, 3.0), "enx")


class TestPermutationLow:
    """Unit tests for the bare ``_permute_ned`` helper."""

    def test_default_enu(self) -> None:
        # NED=(1,2,3) → ENU=(2, 1, -3)
        assert _permute_ned((1.0, 2.0, 3.0), "enu") == (2.0, 1.0, -3.0)

    def test_ned_identity(self) -> None:
        assert _permute_ned((1.0, 2.0, 3.0), "ned") == (1.0, 2.0, 3.0)

    def test_seu(self) -> None:
        # NED=(1,2,3) → SEU=(-N, E, -D) = (-1, 2, -3)
        assert _permute_ned((1.0, 2.0, 3.0), "seu") == (-1.0, 2.0, -3.0)

    def test_all_negations(self) -> None:
        assert _permute_ned((1.0, 2.0, 3.0), "wsu") == (-2.0, -1.0, -3.0)


class TestRotationBasis:
    """Validate :func:`_strike_dip_rake_to_ned` produces a right-handed basis.

    With strike=0 (along north), dip=90 (vertical fault), rake=0 (pure
    strike-slip), the SRF fault-local basis maps into NED as:
      * u1 (along strike)         → (1, 0,  0)   (= north)
      * u2 (orthogonal, in plane) → (0, 0, -1)   (= up; right-handed with u1, u3)
      * u3 (fault normal)         → (0, 1,  0)   (= east)

    Notes on u2: the SRF docs describe u2 as "orthogonal to the strike
    direction but lies in the fault plane" without fixing the sign. The C++
    rconv's rotation formula realises u2 as the *up-dip* direction (so
    {u1, u2, u3} is right-handed with u3 pointing outward), which is what we
    reproduce here verbatim.
    """

    def test_u1_is_strike(self) -> None:
        n, e, d = _strike_dip_rake_to_ned(0, 90, 0, 1, 0, 0)
        assert (n, e, d) == pytest.approx((1.0, 0.0, 0.0), abs=1e-12)

    def test_u2_is_updip_for_vertical_fault(self) -> None:
        n, e, d = _strike_dip_rake_to_ned(0, 90, 0, 0, 1, 0)
        # d = -1 means upwards in NED.
        assert (n, e, d) == pytest.approx((0.0, 0.0, -1.0), abs=1e-12)

    def test_u3_is_normal(self) -> None:
        n, e, d = _strike_dip_rake_to_ned(0, 90, 0, 0, 0, 1)
        assert (n, e, d) == pytest.approx((0.0, 1.0, 0.0), abs=1e-12)

    def test_basis_is_orthonormal(self) -> None:
        # For any strike/dip/rake, {u1, u2, u3} → an orthonormal triple in NED.
        import numpy as np
        for s, d, r in [(0, 90, 0), (150, 90, 0), (122, 40, 90), (45, 30, 60)]:
            vs = [
                _strike_dip_rake_to_ned(s, d, r, *e)
                for e in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            ]
            m = np.array(vs)
            # Should satisfy m @ m.T ≈ I  (i.e. an orthogonal matrix).
            assert np.allclose(m @ m.T, np.eye(3), atol=1e-12)


class TestRotationDocExample:
    """Match the documentation example exactly.

    The SeisSol docs (standard-rupture-format.html) give expected u1/u2/u3 in
    an MCS with ``+axis=ned`` for a vertical strike-slip fault striking 150°
    (the docs show only the resulting vectors, not the source angles; the
    angles below reproduce those vectors exactly under the C++ formula):

        u_1 = {-0.866025403784439,  0.5,                  0}
        u_2 = { 0,                  0,                   -1}
        u_3 = {-0.5,               -0.866025403784439,    0}
    """

    @pytest.fixture
    def mapping(self) -> CoordinateMapping:
        # +axis=ned: output is (north, east, down).
        return CoordinateMapping("+proj=lonlat +datum=WGS84 +units=m +axis=ned")

    def test_u1(self, mapping: CoordinateMapping) -> None:
        v = mapping.to_mcs(150, 90, 0, 1, 0, 0)
        assert v.as_tuple() == pytest.approx(
            (-0.866025403784439, 0.5, 0.0), abs=1e-12,
        )

    def test_u2(self, mapping: CoordinateMapping) -> None:
        v = mapping.to_mcs(150, 90, 0, 0, 1, 0)
        assert v.as_tuple() == pytest.approx((0.0, 0.0, -1.0), abs=1e-12)

    def test_u3(self, mapping: CoordinateMapping) -> None:
        v = mapping.to_mcs(150, 90, 0, 0, 0, 1)
        assert v.as_tuple() == pytest.approx(
            (-0.5, -0.866025403784439, 0.0), abs=1e-12,
        )


class TestToMesh:
    """Geographic → MCS coordinate transforms."""

    def test_projected_tmerc_enu(self) -> None:
        # tmerc with +axis=enu: the result must be (east, north, up).
        # depth_km > 0 means below the surface → z should be negative.
        m = CoordinateMapping(
            "+proj=tmerc +datum=WGS84 +k=0.9996 "
            "+lon_0=-118.5150 +lat_0=34.3440 +units=m +axis=enu"
        )
        v = m.to_mesh(-118.6049, 34.3864, 5.3214)
        # east < 0 (we're west of the origin), north > 0 (we're north of it).
        assert v.x < 0
        assert v.y > 0
        # depth converted to metres and signed for "up": -5321.4 m exactly.
        assert v.z == pytest.approx(-5321.4)

    def test_projected_ned_swaps_xy_and_flips_z(self) -> None:
        # +axis=ned: the (east, north) we got with +enu should swap and
        # the z should flip sign.
        m_enu = CoordinateMapping(
            "+proj=tmerc +datum=WGS84 +k=0.9996 "
            "+lon_0=-118.5150 +lat_0=34.3440 +units=m +axis=enu"
        )
        m_ned = CoordinateMapping(
            "+proj=tmerc +datum=WGS84 +k=0.9996 "
            "+lon_0=-118.5150 +lat_0=34.3440 +units=m +axis=ned"
        )
        v_enu = m_enu.to_mesh(-118.6049, 34.3864, 5.3214)
        v_ned = m_ned.to_mesh(-118.6049, 34.3864, 5.3214)
        assert v_ned.x == pytest.approx(v_enu.y)   # ned.x = north = enu.y
        assert v_ned.y == pytest.approx(v_enu.x)   # ned.y = east  = enu.x
        assert v_ned.z == pytest.approx(-v_enu.z)  # flipped

    def test_geocentric_is_handled(self) -> None:
        # A geocentric MCS is delegated to pyproj entirely. Sanity-check that
        # the output has plausible magnitude (~ Earth radius, 6.3e6 m).
        m = CoordinateMapping("+proj=geocent +datum=WGS84 +units=m +no_defs")
        v = m.to_mesh(-118.6049, 34.3864, 5.3214)
        r = math.sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2)
        assert 6.3e6 < r < 6.4e6
        # Latitude is +34°, so the geocentric Z must be positive.
        assert v.z > 0

    def test_depth_zero_means_surface(self) -> None:
        # With depth=0, z must be 0 regardless of axis (modulo geocentric).
        m = CoordinateMapping(
            "+proj=tmerc +datum=WGS84 +k=0.9996 "
            "+lon_0=0 +lat_0=0 +units=m +axis=enu"
        )
        v = m.to_mesh(0.0, 0.0, 0.0)
        assert v.z == pytest.approx(0.0)


class TestUnitWarnings:
    def test_non_metre_unit_logs_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING, logger="rconv.mapping"):
            # `+proj=lonlat` axes are reported in degrees, which trips the
            # warning.
            CoordinateMapping("+proj=lonlat +datum=WGS84 +units=m +axis=ned")
        warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert any("metre" in r.message for r in warnings)


class TestVec3:
    def test_as_tuple(self) -> None:
        assert Vec3(1.0, 2.0, 3.0).as_tuple() == (1.0, 2.0, 3.0)
