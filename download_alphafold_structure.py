import requests
import os

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

def main():

    uniprot_id = input("Enter UniProt ID: ")
    output_dir = input("Enter output directory (press Enter for current directory): ").strip() or "."
    download_alphafold_structure(uniprot_id, output_dir)\

if __name__ == "__main__":
    main()
