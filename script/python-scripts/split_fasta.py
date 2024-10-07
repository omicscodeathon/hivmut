from Bio import SeqIO
import os
import argparse
from tqdm import tqdm


def split_fasta(input_fasta, output_dir):
    """Splits a multi-FASTA file into individual FASTA files named by their accession IDs.

    Args:
        input_fasta (str): Path to the input multi-FASTA file.
        output_dir (str): Directory to save the individual FASTA files.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Parse the multi-FASTA file and count sequences
    with open(input_fasta, "r") as input_handle:
        fasta_sequences = list(SeqIO.parse(input_handle, "fasta"))

        # Use tqdm to add a progress bar
        for seq_record in tqdm(fasta_sequences, desc="Splitting sequences"):
            # Generate file name using accession ID
            accession_id = seq_record.id
            output_file = os.path.join(output_dir, f"{accession_id}.fasta")

            # Write individual sequence to its own FASTA file
            with open(output_file, "w") as output_handle:
                SeqIO.write(seq_record, output_handle, "fasta")

    print(f"Splitting complete. Individual FASTA files saved in '{output_dir}'.")


def main():
    parser = argparse.ArgumentParser(
        description="Split a multi-FASTA file into individual FASTA files by accession ID."
    )
    parser.add_argument("input_fasta", help="Path to the input multi-FASTA file.")
    parser.add_argument(
        "output_dir", help="Directory where individual FASTA files will be saved."
    )

    args = parser.parse_args()

    # Run the split function
    split_fasta(args.input_fasta, args.output_dir)


if __name__ == "__main__":
    main()
