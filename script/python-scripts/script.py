import argparse
import csv
import os
from tqdm import tqdm
from Bio import Entrez
from Bio import SeqIO


def download_fasta_by_accessions(accession_numbers, output_dir):
    """Downloads FASTA sequences from GenBank based on a list of accession numbers.

    Args:
        accession_numbers (list): A list of accession numbers.
        output_dir (str): The directory to save the downloaded FASTA files.
    """

    Entrez.email = "your_email@example.com"  # Replace with your email address

    for acc_num in tqdm(accession_numbers, desc="Downloading FASTA files"):
        handle = Entrez.efetch(db="nucleotide", id=acc_num, rettype="fasta")
        records = SeqIO.parse(handle, "fasta")
        output_file = f"{output_dir}/{acc_num}.fasta"
        with open(output_file, "w") as f:
            SeqIO.write(records, f, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Download FASTA files from GenBank based on accession numbers."
    )
    parser.add_argument(
        "input_file", help="Path to the input CSV file containing accession numbers."
    )
    parser.add_argument(
        "output_dir", help="Directory to save the downloaded FASTA files."
    )
    args = parser.parse_args()

    # Create the output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Read accession numbers from the CSV file
    with open(args.input_file, "r") as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header row
        accession_numbers = [row[0] for row in reader]

    # Download FASTA files
    download_fasta_by_accessions(accession_numbers, args.output_dir)


if __name__ == "__main__":
    main()
