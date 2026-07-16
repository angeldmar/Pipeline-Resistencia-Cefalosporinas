"""Pruebas unitarias de filter_contigs.py (parte 10: ensamblaje).

Incluye pruebas negativas: FASTA vacio y archivo corrupto (seccion 22 del
diseno del pipeline).
"""

import pytest
from Bio import SeqIO

from filter_contigs import filter_contigs


def write_fasta(path, records: list[tuple[str, int]]) -> None:
    """records: lista de (nombre, longitud) -> escribe secuencias de relleno
    de esa longitud exacta."""
    with open(path, "w") as fasta_file:
        for name, length in records:
            sequence = ("ACGT" * (length // 4 + 1))[:length]
            fasta_file.write(f">{name}\n{sequence}\n")


def test_filter_contigs_keeps_only_long_enough_contigs(tmp_path):
    input_fasta = tmp_path / "contigs.fasta"
    output_fasta = tmp_path / "contigs.filtered.fasta"
    write_fasta(input_fasta, [
        ("contig_long_1", 1200),
        ("contig_short_1", 300),
        ("contig_long_2", 600),
        ("contig_short_2", 100),
        ("contig_exactly_at_threshold", 500),
    ])

    retained_count, total_count = filter_contigs(input_fasta, output_fasta, minimum_length=500)

    assert (retained_count, total_count) == (3, 5)
    retained_names = [record.id for record in SeqIO.parse(output_fasta, "fasta")]
    assert retained_names == ["contig_long_1", "contig_long_2", "contig_exactly_at_threshold"]


def test_filter_contigs_handles_assembly_with_no_contigs_surviving(tmp_path):
    # Caso negativo: un ensamblaje donde NINGUN contig alcanza el minimo no
    # debe hacer fallar el script (results/assemblies/.../contigs.filtered.fasta
    # queda vacio, pero el pipeline sigue -- QUAST luego lo marcara FAIL).
    input_fasta = tmp_path / "contigs.fasta"
    output_fasta = tmp_path / "contigs.filtered.fasta"
    write_fasta(input_fasta, [("only_short_contig", 50)])

    retained_count, total_count = filter_contigs(input_fasta, output_fasta, minimum_length=500)

    assert (retained_count, total_count) == (0, 1)
    assert output_fasta.is_file()
    assert list(SeqIO.parse(output_fasta, "fasta")) == []


def test_filter_contigs_handles_empty_fasta_file(tmp_path):
    # Caso negativo del diseno del pipeline: "FASTQ vacio" -- el equivalente
    # a nivel de ensamblaje es un FASTA de entrada completamente vacio.
    input_fasta = tmp_path / "empty.fasta"
    input_fasta.write_text("")
    output_fasta = tmp_path / "contigs.filtered.fasta"

    retained_count, total_count = filter_contigs(input_fasta, output_fasta, minimum_length=500)

    assert (retained_count, total_count) == (0, 0)


def test_filter_contigs_raises_clear_error_on_corrupted_file(tmp_path):
    # Caso negativo del diseno del pipeline: "archivo corrupto". Un archivo
    # binario/basura que no es texto valido NO debe procesarse en silencio
    # como si tuviera cero contigs validos: debe fallar de forma explicita y
    # temprana (mejor un error claro en este paso que un ensamblaje
    # "vacio pero exitoso" que se descubra recien varios pasos despues).
    corrupted_fasta = tmp_path / "corrupted.fasta"
    corrupted_fasta.write_bytes(bytes(range(256)) * 4)
    output_fasta = tmp_path / "contigs.filtered.fasta"

    with pytest.raises(UnicodeDecodeError):
        filter_contigs(corrupted_fasta, output_fasta, minimum_length=500)


def test_filter_contigs_handles_text_garbage_without_fasta_headers(tmp_path):
    # Variante de "archivo corrupto" que SI es texto valido pero no tiene
    # ningun encabezado FASTA ">": Bio.SeqIO lo interpreta como 0 registros
    # en vez de fallar (comportamiento real verificado, no asumido).
    garbage_fasta = tmp_path / "garbage.fasta"
    garbage_fasta.write_text("esto no es un archivo FASTA valido\nsolo texto suelto\n")
    output_fasta = tmp_path / "contigs.filtered.fasta"

    retained_count, total_count = filter_contigs(garbage_fasta, output_fasta, minimum_length=500)

    assert (retained_count, total_count) == (0, 0)
