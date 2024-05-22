import os.path
import sys
import threading
from Bio import SwissProt,ExPASy
import datetime

def fecth_entrze_id(record):

    entrez_id = ""
    for reference in record.cross_references:

        if ("GeneID" in reference):
            entrez_id = entrez_id + reference[reference.index("GeneID") + 1] + " "
    entrez_id = entrez_id.strip()

    if (len(entrez_id) == 0):
        entrez_id = "None"

    return entrez_id

def fetch_taxonomy_id(record):

    taxonomy_id = ""
    for id in record.taxonomy_id:
        taxonomy_id = taxonomy_id + id + " "

    taxonomy_id = taxonomy_id.strip()

    if (len(taxonomy_id) == 0):
        taxonomy_id = "None"

    return taxonomy_id

def fecth_gene_name(record):

    gene_name_list = record.gene_name
    if(len(gene_name_list)>0 and "Name" in gene_name_list[0]):
        return gene_name_list[0]["Name"]
    else:
        return "None"

def fetch_protein_by_uniprot_id(uniprot_id):
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)

        return record.sequence, datetime.datetime.strptime(record.created[0], "%d-%b-%Y").strftime("%Y-%m-%d"), \
               record.data_class, record.organism, fecth_gene_name(record), fecth_entrze_id(record), fetch_taxonomy_id(record)
               #datetime.datetime.strptime(record.sequence_update[0], "%d-%b-%Y").strftime("%Y-%m-%d"),\
               #datetime.datetime.strptime(record.annotation_update[0], "%d-%b-%Y").strftime("%Y-%m-%d"),\


    except Exception as e:
        print(f"Error fetching protein sequence: {e}")
        return None, None, None, None, None, None, None

def query_single_sequence(): # query a protein sequence with information using UniProt ID

    uniprot_id = input("Enter UniProt ID: ")
    protein_sequence, deposited_date, is_review, organism, gene_name, entrez_id, taxonomy_id = fetch_protein_by_uniprot_id(uniprot_id)
    if protein_sequence:
        print(f"Protein sequence for {uniprot_id}:")
        print(protein_sequence)
        print(f"Deposited date: {deposited_date}")
        #print(f"Sequence update date: {sequence_update}")
        #print(f"Annotation Update date: {annotation_update}")
        print(f"Review/Unreview: {is_review}")
        print(f"Organism: {organism}")
        print(f"Coding-Gene (Name): {gene_name}")
        print(f"Coding-Gene (Entrez ID): {entrez_id}")
        print(f"Taxonomy ID: {taxonomy_id}")

def download_single_sequence(uniprot_id, sequence_file, info_file): # download a protein sequence with information using UniProt ID

    protein_sequence, deposited_date, is_review, organism, gene_name, entrez_id, taxonomy_id = fetch_protein_by_uniprot_id(uniprot_id)

    if(protein_sequence):

        f = open(sequence_file, "w")
        f.write(">" + uniprot_id + "\n" + protein_sequence + "\n")
        f.close()

        f = open(info_file, "w")
        f.write(deposited_date + "\n")
        f.write(is_review + "\n")
        f.write(organism + "\n")
        f.write(gene_name + "\n")
        f.write(entrez_id + "\n")
        f.write(taxonomy_id + "\n")
        f.close()

        print(uniprot_id + " sucessfully download")

    else:
        print(uniprot_id + " failed download")

def download_bulk_sequence(uniprot_id_list, sequence_dir, info_dir): # download multiple protein sequences with information using UniProt ID

    for uniprot_id in uniprot_id_list:

        sequence_file = sequence_dir + "/" + uniprot_id + ".fasta"
        info_file = info_dir + "/" + uniprot_id
        download_single_sequence(uniprot_id, sequence_file, info_file)

def split_name_list(name_list_file, num_sublists):

    f = open(name_list_file, "r")
    text = f.read()
    f.close()

    name_list = text.splitlines()

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

def download_bulk_sequence_multi_thread(uniprot_id_list_file, sequence_dir, info_dir, thread_number): # download protein sequences with information using UniProt ID from a protein list file

    if(os.path.exists(sequence_dir)==False):
        os.makedirs(sequence_dir)

    if (os.path.exists(info_dir) == False):
        os.makedirs(info_dir)

    uniprot_id_list_array = split_name_list(uniprot_id_list_file, thread_number)
    threads = []

    for uniprot_id_list in uniprot_id_list_array:

        thread = threading.Thread(target=download_bulk_sequence, args=(uniprot_id_list, sequence_dir, info_dir,))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

def check_missing_proteins(uniprot_id_list_file, sequence_dir): # check the proteins which failed in fetch_protein_by_uniprot_id()

    f = open(uniprot_id_list_file, "r")
    text = f.read()
    f.close()
    uniprot_id_list = text.splitlines()

    for uniprot_id in uniprot_id_list:
        if(os.path.exists(sequence_dir + "/" + uniprot_id + ".fasta")==False):
            print(uniprot_id)


if __name__ == "__main__":

    #download_bulk_sequence_multi_thread(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    # example
    # python download_sequence_uniprot.py ./all_protein_list ./sequence ./sequence_info 100

    #query_single_sequence()
    # example
    # python download_sequence_uniprot.py

    check_missing_proteins(sys.argv[1], sys.argv[2])
    # example
    # python download_sequence_uniprot.py ./all_protein_list ./sequence/