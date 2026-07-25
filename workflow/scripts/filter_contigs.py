"""Filtrado de contigs cortos en los ensamblajes generados por SPAdes.

SPAdes suele producir, ademas del genoma real, contigs muy cortos que son mas
ruido de ensamblaje (bordes de repeticiones, artefactos de novo) que
secuencia genomica confiable. Este script conserva unicamente los contigs con
al menos una longitud minima configurable, antes de que el ensamblaje pase a
QUAST (evaluacion) o AMRFinderPlus (deteccion de resistencia).
"""

from __future__ import annotations

from pathlib import Path
import argparse

from Bio import SeqIO


def filter_contigs(input_fasta: Path, output_fasta: Path, minimum_length: int = 500) -> tuple[int, int]:
    """Escribe en output_fasta solo los contigs de input_fasta que tengan al
    menos minimum_length bases. Devuelve (contigs_retenidos, contigs_totales)."""
    # "fasta-pearson" en vez del "fasta" a secas: desde Biopython 1.83 el
    # parser "fasta" deprecio el manejo silencioso de texto antes de la primera
    # secuencia y en una version futura lanzara ValueError ante el. Para un
    # ensamblaje real de SPAdes (cabeceras ">NODE_...") el resultado es
    # identico, y "fasta-pearson" preserva el manejo gracioso del que dependen
    # las pruebas negativas: un archivo de texto sin ninguna cabecera ">" se
    # interpreta como cero contigs (que el pipeline reporta como ensamblaje
    # vacio) en lugar de hacer fallar el script. Un archivo binario/no-texto
    # sigue fallando temprano con UnicodeDecodeError, como debe.
    all_contigs = list(SeqIO.parse(input_fasta, "fasta-pearson"))
    retained_contigs = [contig for contig in all_contigs if len(contig.seq) >= minimum_length]

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(retained_contigs, output_fasta, "fasta")

    return len(retained_contigs), len(all_contigs)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filtrar contigs cortos de un ensamblaje FASTA."
    )
    parser.add_argument("input_fasta", type=Path, help="Ensamblaje FASTA sin filtrar (ej. contigs.fasta de SPAdes)")
    parser.add_argument("output_fasta", type=Path, help="Ruta de salida para el FASTA filtrado")
    parser.add_argument(
        "--minimum-length", type=int, default=500,
        help="Longitud minima (pb) para conservar un contig (config: assembly.minimum_contig_length)",
    )
    args = parser.parse_args()

    retained_count, total_count = filter_contigs(args.input_fasta, args.output_fasta, args.minimum_length)
    discarded_count = total_count - retained_count

    if retained_count == 0:
        print(
            f"ADVERTENCIA: {args.input_fasta.name} se quedo sin contigs tras el filtrado "
            f"(0/{total_count} >= {args.minimum_length} pb). Revisar el ensamblaje."
        )
    else:
        print(
            f"{args.input_fasta.name}: {retained_count}/{total_count} contigs retenidos "
            f"(>= {args.minimum_length} pb); {discarded_count} descartados."
        )


if __name__ == "__main__":
    main()
