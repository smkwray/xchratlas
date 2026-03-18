from xatlas.region import classify_interval_region, classify_x_region


def test_classify_x_region_boundaries():
    assert classify_x_region(60001) == "PAR1"
    assert classify_x_region(2699520) == "PAR1"
    assert classify_x_region(2699521) == "nonPAR"
    assert classify_x_region(154931043) == "nonPAR"
    assert classify_x_region(154931044) == "PAR2"
    assert classify_x_region(155260560) == "PAR2"


def test_classify_interval_region():
    assert classify_interval_region(70000, 80000) == "PAR1"
    assert classify_interval_region(3_000_000, 4_000_000) == "nonPAR"
    assert classify_interval_region(155_000_000, 155_100_000) == "PAR2"
    assert classify_interval_region(2_699_000, 2_700_500) == "boundary_crossing"
