import sys

import requests
import os
import threading

def download_alphafold_structure(uniprot_id, output_dir="."):
    try:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id.upper()}-F1-model_v4.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            output_file = os.path.join(output_dir, f"{uniprot_id}.pdb")
            with open(output_file, "wb") as f:
                f.write(response.content)
            print(f"AlphaFold2 structure for {uniprot_id} downloaded successfully.")
            print(f"Saved to: {output_file}")
        else:
            print(f"Failed to retrieve AlphaFold2 structure for {uniprot_id}.")
    except Exception as e:
        print(f"Error downloading AlphaFold2 structure: {e}")

def download_single_structure_by_input_uniprot_id(): # download a protein structure from AF2 database by inputting UniProt ID and save dir

    uniprot_id = input("Enter UniProt ID: ")
    output_dir = input("Enter output directory (press Enter for current directory): ").strip() or "."
    download_alphafold_structure(uniprot_id, output_dir)

def download_bulk_structure(uniprot_id_list, structure_dir):

    for uniprot_id in uniprot_id_list:
        download_alphafold_structure(uniprot_id, structure_dir)

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

def download_bulk_structure_multi_thread(uniprot_id_list_file, structure_dir, thread_number): # download protein sequences with information using UniProt ID from a protein list file

    if (os.path.exists(structure_dir) == False):
        os.makedirs(structure_dir)

    uniprot_id_list_array = split_name_list(uniprot_id_list_file, thread_number)
    threads = []

    for uniprot_id_list in uniprot_id_list_array:

        thread = threading.Thread(target=download_bulk_structure, args=(uniprot_id_list, structure_dir,))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()


def check_missing_proteins(uniprot_id_list_file, structure_dir): # check the proteins which failed in download_alphafold_structure()

    f = open(uniprot_id_list_file, "r")
    text = f.read()
    f.close()
    uniprot_id_list = text.splitlines()

    for uniprot_id in uniprot_id_list:
        if(os.path.exists(structure_dir + "/" + uniprot_id + ".pdb")==False):
            print(uniprot_id)


if __name__ == "__main__":

    #download_single_structure_by_input_uniprot_id()
    #example: python download_alphafold_structure.py

    download_bulk_structure_multi_thread(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    #example: python download_alphafold_structure.py test_name_list ./AF2_structure/ 100

    #check_missing_proteins(sys.argv[1], sys.argv[2])


