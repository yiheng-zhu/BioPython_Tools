import sys
import threading
from Bio import SwissProt,ExPASy
import os

def fetch_protein_by_uniprot_id(uniprot_id, reference_file):
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)

        for i in range(len(record.cross_references)):
            reference = record.cross_references[i]
            if('RefSeq' in reference):
                f = open(reference_file, "w")
                for i in range(1, len(reference)):
                    f.write(reference[i] + "\n")
                f.close()
                break

    except Exception as e:
        print(f"Error fetching protein sequence: {e}")

def download_bulk_sequence(uniprot_id_list, reference_dir):

    for uniprot_id in uniprot_id_list:

        reference_file = reference_dir + "/" + uniprot_id
        fetch_protein_by_uniprot_id(uniprot_id, reference_file)

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

def download_bulk_sequence_multi_thread(uniprot_id_list_file, reference_dir, thread_number): # download protein sequences with information using UniProt ID from a protein list file


    if (os.path.exists(reference_dir) == False):
        os.makedirs(reference_dir)

    uniprot_id_list_array = split_name_list(uniprot_id_list_file, thread_number)
    threads = []

    for uniprot_id_list in uniprot_id_list_array:

        thread = threading.Thread(target=download_bulk_sequence, args=(uniprot_id_list, reference_dir,))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

def check_missing_proteins(uniprot_id_list_file, reference_dir): # check the proteins which failed in fetch_protein_by_uniprot_id()

    f = open(uniprot_id_list_file, "r")
    text = f.read()
    f.close()
    uniprot_id_list = text.splitlines()

    for uniprot_id in uniprot_id_list:
        if(os.path.exists(reference_dir + "/" + uniprot_id )==False):
            print(uniprot_id)

if __name__ == "__main__":

    download_bulk_sequence_multi_thread(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    #check_missing_proteins(sys.argv[1], sys.argv[2])