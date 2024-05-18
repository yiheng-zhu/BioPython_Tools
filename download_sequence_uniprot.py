from Bio import SwissProt,ExPASy

def fetch_protein_sequence(uniprot_id):
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)
        return record.sequence

    except Exception as e:
        print(f"Error fetching protein sequence: {e}")
        return None

def main():
    uniprot_id = input("Enter UniProt ID: ")
    protein_sequence = fetch_protein_sequence(uniprot_id)
    if protein_sequence:
        print(f"Protein sequence for {uniprot_id}:")
        print(protein_sequence)


if __name__ == "__main__":
    main()