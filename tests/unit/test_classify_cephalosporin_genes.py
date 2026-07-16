"""Pruebas unitarias de classify_cephalosporin_genes.py (parte 15).

El caso critico de esta suite: confirmar que NO se asume que todo TEM/SHV
detectado es una BLEE solo por el prefijo del simbolo (seccion 13.3 del
diseno del pipeline).
"""

import yaml

from classify_cephalosporin_genes import classify_beta_lactamase, load_resistance_targets


def load_real_resistance_targets(repo_root):
    return load_resistance_targets(repo_root / "config" / "resistance_targets.yaml")


def test_ctx_m_is_always_esbl_regardless_of_subclass(repo_root):
    targets = load_real_resistance_targets(repo_root)
    assert classify_beta_lactamase("blaCTX-M-15", "CEPHALOSPORIN", targets) == "ESBL"


def test_tem_with_extended_spectrum_subclass_is_esbl(repo_root):
    targets = load_real_resistance_targets(repo_root)
    assert classify_beta_lactamase("blaTEM-52", "EXTENDED-SPECTRUM BETA-LACTAM", targets) == "ESBL"


def test_tem_without_extended_spectrum_subclass_is_not_esbl(repo_root):
    # blaTEM-1 es la penicilinasa clasica de espectro estrecho, no una BLEE.
    # Este es el caso central que el documento pide no asumir mal.
    targets = load_real_resistance_targets(repo_root)
    result = classify_beta_lactamase("blaTEM-1", "BETA-LACTAM", targets)
    assert result != "ESBL"
    assert result == "Other"


def test_shv_without_extended_spectrum_subclass_is_not_esbl(repo_root):
    targets = load_real_resistance_targets(repo_root)
    result = classify_beta_lactamase("blaSHV-1", "BETA-LACTAM", targets)
    assert result != "ESBL"


def test_cmy_is_ampc(repo_root):
    targets = load_real_resistance_targets(repo_root)
    assert classify_beta_lactamase("blaCMY-2", "CEPHALOSPORIN", targets) == "AmpC"


def test_ndm_is_carbapenemase(repo_root):
    targets = load_real_resistance_targets(repo_root)
    assert classify_beta_lactamase("blaNDM-1", "CARBAPENEM", targets) == "Carbapenemase"


def test_unrecognized_gene_is_other(repo_root):
    targets = load_real_resistance_targets(repo_root)
    assert classify_beta_lactamase("blaXYZ-99", "BETA-LACTAM", targets) == "Other"


def test_load_resistance_targets_reads_real_config_file(repo_root):
    targets = load_real_resistance_targets(repo_root)
    assert "CTX-M" in targets["extended_spectrum_beta_lactamases"]
    assert "CMY" in targets["ampc"]
    assert "NDM" in targets["carbapenemases"]
