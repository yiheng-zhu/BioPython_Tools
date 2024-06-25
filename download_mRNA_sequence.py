import os
import sys
import threading
from Bio import Entrez, SeqIO


def fetch_mrna_id(entrez_id):

    try:

        Entrez.email = "your@email.com"  # Provide your email to NCBI

        handle = Entrez.efetch(db="gene", id=entrez_id, retmode="xml")
        record_list = Entrez.read(handle)
        handle.close()
        for record in record_list:
            if ('Entrezgene_locus' in record):
                sub_record = record['Entrezgene_locus']
                for sub_sub_record in sub_record:
                    if('Gene-commentary_products' in sub_sub_record):
                        element_list = sub_sub_record['Gene-commentary_products']
                        for element in element_list:
                            product_type = element['Gene-commentary_type'].attributes["value"]
                            product_id = element['Gene-commentary_accession']
                            if (product_type == "mRNA"):
                                return product_id

        return ""

    except Exception as e:
        print(f"Error mapping Entrez ID {entrez_id} to GI number with regions: {e}")
        return ""

def fetch_mrna_sequence(entrez_id):
    # Provide your email address to NCBI Entrez
    try:
        Entrez.email = "your@email.com"

        mrna_id = fetch_mrna_id(entrez_id)
        print(mrna_id)
        if(len(mrna_id)>0):

            # Fetch the gene record from GenBank
            handle = Entrez.efetch(db="nucleotide", id=mrna_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            sequence = record.seq
            return sequence
        else:
            return ""

    except Exception as e:
        print(f"Error mapping Entrez ID {entrez_id} to GI number with regions: {e}")
        return ""


def fetch_mrna_sequence_by_name(reference_file, sequence_file):


    try:

        f = open(reference_file, "r")
        f.readline()
        line = f.readline()
        name = line.split()[0].strip(".")

        Entrez.email = "your@email.com"

        handle = Entrez.efetch(db="nucleotide", id = name, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        sequence = record.seq
        if(len(sequence)>=0):
            f = open(sequence_file, "w")
            f.write(">" + name + "\n" + str(sequence) + "\n")
            f.close()

    except Exception as e:
        print(f"Error mapping Entrez ID {reference_file} to GI number with regions: {e}")

def download_bulk_sequence(reference_dir, mrna_sequence_dir, name_list):

    for name in name_list:
        reference_file = reference_dir + "/" + name
        mrna_sequence_file = mrna_sequence_dir + "/" + name

        if(os.path.exists(mrna_sequence_file)==False):
            fetch_mrna_sequence_by_name(reference_file, mrna_sequence_file)

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

def download_bulk_sequence_multi_thread(reference_dir, mrna_sequence_dir, thread_number): # download protein sequences with information using UniProt ID from a protein list file

    if(os.path.exists(mrna_sequence_dir)==False):
        os.makedirs(mrna_sequence_dir)

    all_name_list = os.listdir(reference_dir)

    name_list_array = split_name_list(all_name_list, thread_number)
    threads = []

    for name_list in name_list_array:

        thread = threading.Thread(target=download_bulk_sequence, args=(reference_dir, mrna_sequence_dir, name_list, ))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

if __name__ == "__main__":

    download_bulk_sequence_multi_thread(sys.argv[1], sys.argv[2], int(sys.argv[3]))

