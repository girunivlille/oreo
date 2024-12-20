import numpy as np
import matplotlib.pyplot as plt
import argparse
import math
import sys

# Lit les trois fichiers de benchmark (dans l'ordre untouched, shuffled et sorted)
# Produit un graphique par mesure (temps et mémoire) pour les 3 fichiers

def graph_from_memtimes(log_file, pie_graph_file):
    with open(log_file, 'r') as log:
        lines = log.readlines()
    unmapped = int(lines[0].strip().split(':')[1].strip())
    mapped = int(lines[1].strip().split(':')[1].strip())
    still_unmapped = int(lines[2].strip().split(':')[1].strip())
    all_reads = int(lines[3].strip().split(':')[1].strip())

    in_contigs = all_reads - unmapped

    pie_data = [in_contigs, mapped, still_unmapped]

    plt.pie(pie_data, labels = ["in_contigs", "linked", "unmapped"], autopct=lambda pct: str(math.floor(all_reads*pct/100)))
    plt.title("Répartition des reads dans leur catégorie de mapping")
    plt.savefig(pie_graph_file)


log_file = sys.argv[1]
pie_graph_file = sys.argv[2]

graph_from_memtimes(log_file, pie_graph_file)
