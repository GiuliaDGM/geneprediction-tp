# IMPORTS

import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     f"{sys.argv[0]} -h")
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path,
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


###########################################################
############# Lecture du fichier de séquences #############
###########################################################


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    with fasta_file.open('r') as file:
        lines = file.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence.upper()


###################################################################################################
######### Recherche des codons d’initiation, stop et de la séquence de Shine-Dalgarno #############
###################################################################################################


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    match = start_regex.search(sequence, start, stop)
    return match.start(0) if match else None


def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    for match in stop_regex.finditer(sequence, start):
        if (match.start(0) - start) % 3 == 0:
            return match.start(0)
    return None


def has_shine_dalgarno(shine_regex: Pattern,
                       sequence: str,
                       start: int,
                       max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno 
                                        to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    # Calculer la plage de recherche
    search_start = start - max_shine_dalgarno_distance
    search_end = start - 6

    # Vérifier si la position de début de recherche est valide
    if search_start < 0:
        return False

    # Rechercher le motif Shine-Dalgarno dans l'intervalle
    if shine_regex.search(sequence, search_start, search_end):
        return True
    return False

def predict_genes(sequence: str, start_regex: Pattern, stop_regex: Pattern, shine_regex: Pattern,
                  min_gene_len: int, max_shine_dalgarno_distance: int, min_gap: int) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno 
                                        to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    genes = []
    position_courante = 0
    while len(sequence) - position_courante >= min_gap:
        start = find_start(start_regex, sequence, position_courante, len(sequence))
        if start is None:
            break
        stop = find_stop(stop_regex, sequence, start)
        if stop and (stop - start >= min_gene_len):
            if has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
                genes.append([start + 1, stop + 3])  # Positions de début et fin en format 1-based
                position_courante = stop + 3 + min_gap
            else:
                position_courante += 1
        else:
            position_courante += 1
    return genes


def write_genes_pos(predicted_genes_file: Path, probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit(f"Error cannot open {predicted_genes_file}")


def write_genes(fasta_file: Path, sequence: str, probable_genes: List[List[int]], sequence_rc: str,
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position 
                                of each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt", encoding="utf-8") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(f">gene_{i+1}{os.linesep}{textwrap.fill(sequence[gene_pos[0]-1:gene_pos[1]])}{os.linesep}")
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(f">gene_{i+1+j}{os.linesep}{textwrap.fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])}{os.linesep}")
    except IOError:
        sys.exit(f"Error cannot open {fasta_file}")


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


#==============================================================
# Main program
#==============================================================


def main() -> None: # pragma: no cover
    """
    Main program function
    """
    # Initialisation des regex pour les codons start, stop et Shine-Dalgarno
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')

    # Arguments
    args = get_arguments()

    # Chargement de la séquence du génome
    sequence = read_fasta(args.genome_file)

    # Prédiction des gènes dans le sens 5' vers 3'
    probable_genes = predict_genes(
        sequence,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap
    )

    # Génération du complément inverse de la séquence pour le sens 3' vers 5'
    sequence_rc = reverse_complement(sequence)

    # Prédiction des gènes dans le sens 3' vers 5' sur la séquence complémentaire
    probable_genes_comp = predict_genes(
        sequence_rc,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap
    )

    # Correction des positions des gènes trouvés dans le sens 3' vers 5'
    corrected_genes_comp = [[len(sequence) - end + 1, len(sequence) - start + 1] for start,
                            end in probable_genes_comp]

    # Écriture des positions des gènes dans le fichier de sortie
    write_genes_pos(args.predicted_genes_file, probable_genes + corrected_genes_comp)

    # Écriture des séquences des gènes dans le fichier FASTA de sortie
    write_genes(args.fasta_file, sequence, probable_genes, sequence_rc, corrected_genes_comp)


if __name__ == '__main__':
    main()
