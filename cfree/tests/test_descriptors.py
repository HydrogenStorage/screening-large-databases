from math import isclose

from cfree.descriptors import compute_wth2, saturate_molecule, count_h2_difference


def test_hydrogenate():
    assert saturate_molecule('C=C') == 'CC'
    assert saturate_molecule('c1cnccc1') == 'C1CNCCC1'


def test_wth2():
    assert isclose(compute_wth2('C=C=C'), 9.14, abs_tol=1e-2)
    assert isclose(compute_wth2('C=C=C', 'C=CC'), 4.79, abs_tol=1e-2)


def test_count():
    assert count_h2_difference('C=C') == 1
    assert count_h2_difference('CC') == 0
