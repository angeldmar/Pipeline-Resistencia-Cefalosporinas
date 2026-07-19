"""Pruebas unitarias de parse_abricate.py (parte 25: segundo motor de AMR)."""

import pandas as pd
import pytest

from parse_abricate import derive_gene_family, normalize_abricate_table


def build_raw_abricate_row(gene="blaCTX-M-15", identity=99.78, coverage=100.0, database="card") -> dict:
    return {
        "FILE": "contigs.fasta",
        "SEQUENCE": "contig_1",
        "START": 100,
        "END": 1000,
        "STRAND": "+",
        "GENE": gene,
        "COVERAGE": "1-900/900",
        "COVERAGE_MAP": "===============",
        "GAPS": "0/0",
        "%COVERAGE": coverage,
        "%IDENTITY": identity,
        "DATABASE": database,
        "ACCESSION": "ARO:3002225",
        "PRODUCT": "CTX-M-15 beta-lactamase",
        "RESISTANCE": "cephalosporin",
    }


def test_derive_gene_family_strips_allele_suffix():
    assert derive_gene_family("blaCTX-M-15") == "blaCTX-M"
    assert derive_gene_family("blaTEM-1") == "blaTEM"
    assert derive_gene_family("blaKPC") == "blaKPC"  # sin sufijo numerico, queda igual


def test_derive_gene_family_strips_chained_suffixes():
    # ResFinder agrega su propio sufijo de desambiguacion interno ademas del
    # alelo (formato real observado corriendo ABricate/ResFinder real sobre
    # una muestra con blaCMY-2): "blaCMY-2_1" debe reducirse a "blaCMY", no
    # quedar a medias en "blaCMY-2" -- que rompería la coincidencia de
    # familia contra el "blaCMY" que entrega AMRFinderPlus para el mismo gen.
    assert derive_gene_family("blaCMY-2_1") == "blaCMY"


def test_normalize_abricate_table_computes_confidence_threshold():
    raw_table = pd.DataFrame([
        build_raw_abricate_row(gene="blaCTX-M-15", identity=99.78, coverage=100.0),
        build_raw_abricate_row(gene="blaTEM-1", identity=85.0, coverage=60.0, database="resfinder"),
    ])

    normalized_table = normalize_abricate_table(
        "EC001", raw_table, minimum_identity=90, minimum_gene_coverage=80
    )

    assert len(normalized_table) == 2
    ctx_m_row = normalized_table.loc[normalized_table["gene_symbol"] == "blaCTX-M-15"].iloc[0]
    tem_row = normalized_table.loc[normalized_table["gene_symbol"] == "blaTEM-1"].iloc[0]

    assert ctx_m_row["gene_family"] == "blaCTX-M"
    assert ctx_m_row["meets_identity_coverage_threshold"] == True  # noqa: E712 (numpy bool)
    assert tem_row["meets_identity_coverage_threshold"] == False  # noqa: E712 (60% cobertura < 80%)


def test_normalize_abricate_table_handles_empty_input():
    empty_raw_table = pd.DataFrame(columns=[
        "FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE", "COVERAGE_MAP",
        "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE",
    ])

    normalized_table = normalize_abricate_table(
        "EC_NO_HITS", empty_raw_table, minimum_identity=90, minimum_gene_coverage=80
    )

    assert normalized_table.empty
