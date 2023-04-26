from cfree.compute import compute_storage_energy


def test_compute():
    eng = compute_storage_energy('C=C')
    assert eng > 0
