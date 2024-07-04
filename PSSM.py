import subprocess
import numpy as np
import sys
from decimal import Decimal
import math
import os
import threading

blast_file = "/data/yihengzhu/toolbars/sequence_homology_tools/PSSM/ncbi-blast-2.15.0+/bin/"
db_file = "/data/yihengzhu/toolbars/sequence_homology_tools/PSSM/uniref50_database/uniref50"
cut_off = 1024

def run_psi_blast(query_file, output_file, iterations = 3, evalue = 0.001):
    cmd = [
        blast_file + "/psiblast",
        "-query", query_file,
        "-db", db_file,
        "-num_iterations", str(iterations),
        "-evalue", str(evalue),
        "-out_ascii_pssm", output_file,
        "-outfmt", "0"  # output format 0 for the traditional BLAST output
    ]

    # Run the PSI-BLAST command
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running PSI-BLAST: {e}")


def normalize_pssm(origin_pssm_file, log_pssm_file):

    f = open(log_pssm_file, "w")
    pssm_matrix = read_pssm(origin_pssm_file)
    for row in pssm_matrix:
        line = ""
        for element in row:
            value = 1.0/(1.0 + math.exp(-element))
            line = line + str(Decimal(value).quantize(Decimal("0.000000"))) + " "
        line = line.strip()
        f.write(line + "\n")
    f.close()



def read_pssm(pssm_file):

    with open(pssm_file) as f:
        lines = f.readlines()

    # Identify the start of the PSSM matrix
    header_end = 3
    data_start = 3
    while len(lines[header_end].strip())>0:
        header_end += 1

    data_lines = lines[data_start:header_end + 1]

    pssm_matrix = []
    for line in data_lines:
        if not line.strip():
            continue
        parts = line.split()
        scores = [float(x) for x in parts[2:22]]
        pssm_matrix.append(scores)

    return pssm_matrix

def run_single_pssm(query_file, origin_pssm_file, log_pssm_file):

    if(os.path.exists(log_pssm_file)):
        return

    run_psi_blast(query_file, origin_pssm_file)
    if(os.path.exists(origin_pssm_file)):
        normalize_pssm(origin_pssm_file, log_pssm_file)

def run_thread_pssm(sequence_dir, origin_pssm_dir, log_pssm_dir, name_list):

    for name in name_list:

        query_file = sequence_dir + "/" + name
        origin_pssm_file = origin_pssm_dir + "/" + name.split(".")[0] + ".pssm"
        log_pssm_file = log_pssm_dir + "/" + name.split(".")[0] + ".pssm"
        run_single_pssm(query_file, origin_pssm_file, log_pssm_file)

def split_sequence(sequence_file, sequence_dir):

    f = open(sequence_file, "r")
    text = f.read()
    f.close()

    line_set = text.splitlines()
    line_set = line_set[130000:]

    for line in line_set:
        line = line.strip()
        if (line.startswith(">")):
            name = line[1:]
        else:

            if(len(line)>cut_off):
                line = line[0:cut_off]

            f = open(sequence_dir + "/" + name + ".fasta", "w")
            f.write(">" + name + "\n" + line + "\n")
            f.close()

def split_name_list(name_list, num_sublists):

    sublist_length = len(name_list) // num_sublists
    remainder = len(name_list) % num_sublists

    # Create the sublists
    sublists = []
    start = 0
    for i in range(num_sublists):
        if i < remainder:
            end = start + sublist_length + 1
        else:
            end = start + sublist_length
        sublists.append(name_list[start:end])
        start = end

    return sublists

def create_dir(temp_dir):

    if (os.path.exists(temp_dir) == False):
        os.makedirs(temp_dir)


def create_pssm(sequence_file, sequence_dir, origin_pssm_dir, log_pssm_dir, thread_number): # download protein sequences with information using UniProt ID from a protein list file

    create_dir(sequence_dir)
    create_dir(origin_pssm_dir)
    create_dir(log_pssm_dir)

    split_sequence(sequence_file, sequence_dir)

    all_name_list = os.listdir(sequence_dir)

    name_list_array = split_name_list(all_name_list, thread_number)
    threads = []

    for name_list in name_list_array:

        thread = threading.Thread(target=run_thread_pssm, args=(sequence_dir, origin_pssm_dir, log_pssm_dir, name_list, ))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

if __name__ == "__main__":

    create_pssm(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))