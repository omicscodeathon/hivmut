from Bio import SeqIO
import os
import argparse
from tqdm import tqdm


def remove_gaps(seq):
    """Remove gaps (underscores) from a sequence."""
    return seq.replace("-", "")


def process_fasta_dir(input_dir, output_fasta):
    """Processes a directory of FASTA files, removes gaps, and combines the sequences.

    Args:
        input_dir (str): Path to the directory containing FASTA files.
        output_fasta (str): Path to the output FASTA file to save cleaned sequences.
    """
    cleaned_sequences = []
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]

    with tqdm(
        total=len(fasta_files), desc="Processing FASTA files", unit="file"
    ) as pbar:
        for fasta_file in fasta_files:
            fasta_path = os.path.join(input_dir, fasta_file)
            with open(fasta_path, "r") as handle:
                fasta_sequences = SeqIO.parse(handle, "fasta")

                for seq_record in fasta_sequences:
                    # Remove gaps from the sequence
                    seq_record.seq = remove_gaps(str(seq_record.seq))
                    cleaned_sequences.append(seq_record)

            pbar.update(1)

    # Write the cleaned sequences to the output file
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(cleaned_sequences, output_handle, "fasta")

    print(
        f"Processing complete. {len(cleaned_sequences)} sequences saved to {output_fasta}."
    )


def main():
    parser = argparse.ArgumentParser(
        description="Process a directory of FASTA files, remove gaps, and save cleaned sequences."
    )
    parser.add_argument(
        "input_dir", help="Path to the input directory containing FASTA files."
    )
    parser.add_argument(
        "output_fasta", help="Path to the output FASTA file to save cleaned sequences."
    )

    args = parser.parse_args()

    # Run the processing function
    process_fasta_dir(args.input_dir, args.output_fasta)


if __name__ == "__main__":
    main()
