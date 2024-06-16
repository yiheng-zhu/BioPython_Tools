import math
import sys

from Bio import Entrez
from Bio import SeqIO
import os
import threading

def complement_sequence(dna_sequence):

    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'K':'M', 'W':'W', 'M':'K', 'S':'S', 'Y':'R', 'R':'Y'}
    complement = ""
    for i in range(len(dna_sequence)-1, 0, -1):
        complement = complement + complement_map[dna_sequence[i]]

    return complement

def map_entrez_to_gi(entrez_id):
    try:
        # Set up Entrez
        Entrez.email = "your@email.com"  # Provide your email to NCBI

        # Fetch the record for the given Entrez ID
        handle = Entrez.efetch(db="gene", id=entrez_id, retmode="xml")
        record_list = Entrez.read(handle)
        handle.close()
        gi_number = ""
        strand = ""
        region_start = -1
        region_end = -1

        for record in record_list:
            if('Entrezgene_locus' in record):
                sub_record = record['Entrezgene_locus']
                for sub_sub_record in sub_record:
                    if ('Gene-commentary_accession' in sub_sub_record):
                        gi_number = sub_sub_record['Gene-commentary_accession']
                        if('Gene-commentary_seqs' in sub_sub_record):
                            for locate in sub_sub_record['Gene-commentary_seqs']:
                                if('Seq-loc_int' in locate and 'Seq-interval' in locate['Seq-loc_int'] and 'Seq-interval_from' in locate['Seq-loc_int']['Seq-interval'] and 'Seq-interval_to' in locate['Seq-loc_int']['Seq-interval']):
                                    region_start = locate['Seq-loc_int']['Seq-interval']['Seq-interval_from']
                                    region_end = locate['Seq-loc_int']['Seq-interval']['Seq-interval_to']
                                    try:
                                        strand = locate['Seq-loc_int']['Seq-interval']['Seq-interval_strand'][
                                            'Na-strand'].attributes['value']
                                    except Exception as e:
                                        strand = "plus"

                                    if(strand=="minus"):
                                        return gi_number, int(region_start), int(region_end) + 1, strand
                                    else:
                                        return gi_number, int(region_start) + 1, int(region_end) + 1, strand
                                else:
                                    try:
                                        start_list = []
                                        end_list = []
                                        strand_list = []

                                        all_str_line = str(locate)

                                        while("Seq-interval_from" in all_str_line):

                                            start = all_str_line.find("Seq-interval_from")
                                            end = all_str_line.find("Seq-interval_id")
                                            str_line = all_str_line[start:end]
                                            values = str_line.split(",")


                                            start_list.append(int(values[0].strip("'").split(":")[1].strip().strip("'")))
                                            end_list.append(int(values[1].strip("'").split(":")[1].strip().strip("'")))
                                            pos = str_line.find("attributes={'value': ")
                                            pos_line = str_line[pos:]
                                            strand = pos_line.split(":")[1].strip().split("}")[0].strip().strip("'")
                                            all_str_line = all_str_line[end+len("Seq-interval_id"):]
                                            strand_list.append(strand)

                                        all_pos = start_list + end_list
                                        max_value = max(all_pos)
                                        min_value = min(all_pos)

                                        if (strand_list[0] == "minus"):
                                            return gi_number, int(min_value) , int(max_value) + 1, strand_list[0]
                                        else:
                                            return gi_number, int(min_value) + 1, int(max_value) + 1, strand_list[0]
                                    except Exception as e:
                                        print(e)

        return gi_number, region_start, region_end, strand

    except Exception as e:
        print(f"Error mapping Entrez ID {entrez_id} to GI number with regions: {e}")
        return None, None, None, None

def extract_gene_sequence(entrez_id, gi_number, start, end, strand):
    try:
        # Set up Entrez
        Entrez.email = "yihzhu@njau.edu.cn"  # Provide your email to NCBI

        # Fetch the sequence for the given Entrez ID
        handle = Entrez.efetch(db="Nucleotide", id=gi_number, rettype="fasta", retmode="text", seq_start = start, seq_stop = end)

        # Read the sequence
        record = SeqIO.read(handle, "fasta")

        # Close the handle
        handle.close()

        if(strand=="minus"):
            return complement_sequence(record.seq)
        else:
            return str(record.seq)

    except Exception as e:
        print(f"Error extracting sequence for Entrez ID {entrez_id}: {e}")
        return None

def download_single_sequence(entrez_id, sequence_file):

    gi_number, region_start, region_end, strand = map_entrez_to_gi(entrez_id)

    sequence = ""
    if(region_start!=-1):
        sequence = extract_gene_sequence(entrez_id, gi_number, region_start, region_end, strand)
    else:
        print(entrez_id)

    if(sequence):
        f = open(sequence_file, "w")
        f.write(">" + entrez_id + "\n" + sequence + "\n")
        f.close()

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

def download_bulk_sequence(entrez_id_list, sequence_dir):

    for entrez_id in entrez_id_list:
        if(os.path.exists(sequence_dir + "/" + entrez_id + ".fasta")==False):
            download_single_sequence(entrez_id, sequence_dir + "/" + entrez_id + ".fasta")


def download_bulk_sequence_multi_thread(entrez_id_list_file, sequence_dir, thread_number): # download protein sequences with information using UniProt ID from a protein list file

    if(os.path.exists(sequence_dir)==False):
        os.makedirs(sequence_dir)

    entrez_id_list_array = split_name_list(entrez_id_list_file, thread_number)
    threads = []

    for entrez_id_list in entrez_id_list_array:

        thread = threading.Thread(target=download_bulk_sequence, args=(entrez_id_list, sequence_dir, ))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

def check_missing_genes(gene_id_list_file, sequence_dir): # check the proteins which failed in fetch_protein_by_uniprot_id()

    f = open(gene_id_list_file, "r")
    text = f.read()
    f.close()
    gene_id_list = text.splitlines()

    for gene_id in gene_id_list:
        if(os.path.exists(sequence_dir + "/" + gene_id + ".fasta")==False):
            print(gene_id)



if __name__ == "__main__":

    download_bulk_sequence_multi_thread(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    # example: python download_gene_sequence.py ./all_gene_list ./gene_sequence/ 5

    #check_missing_genes(sys.argv[1], sys.argv[2])
    # example: python ./download_gene_sequence.py ./all_gene_list ./gene_sequence/

    #download_single_sequence(sys.argv[1], sys.argv[2])


