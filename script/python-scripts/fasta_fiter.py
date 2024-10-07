from Bio import SeqIO
import os
import argparse
from tqdm import tqdm


def filter_fasta_dir(input_dir, output_fasta, patterns):
    """Filters sequences from FASTA files in a directory based on keywords in the header or sequence.

    Args:
        input_dir (str): Path to the directory containing FASTA files.
        output_fasta (str): Path to the output FASTA file to save filtered sequences.
        patterns (list of str): List of keywords to search for in the header or sequence.
    """
    matched_sequences = []
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]

    with tqdm(
        total=len(fasta_files), desc="Processing FASTA files", unit="file"
    ) as pbar:
        for fasta_file in fasta_files:
            fasta_path = os.path.join(input_dir, fasta_file)
            with open(fasta_path, "r") as handle:
                fasta_sequences = SeqIO.parse(handle, "fasta")

                for seq_record in fasta_sequences:
                    # Check if all patterns are in the sequence description or sequence
                    if all(
                        pattern.lower() in seq_record.description.lower()
                        for pattern in patterns
                    ):
                        matched_sequences.append(seq_record)

            pbar.update(1)

    # Write the matched sequences to the output file
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(matched_sequences, output_handle, "fasta")

    print(
        f"Filtering complete. {len(matched_sequences)} sequences saved to {output_fasta}."
    )


def main():
    parser = argparse.ArgumentParser(
        description="Filter a directory of FASTA files for sequences containing specific patterns."
    )
    parser.add_argument(
        "input_dir", help="Path to the input directory containing FASTA files."
    )
    parser.add_argument(
        "output_fasta", help="Path to the output FASTA file to save filtered sequences."
    )
    parser.add_argument(
        "patterns", nargs="+", help="Keywords to search for (e.g., HIV1 pol)."
    )

    args = parser.parse_args()

    # Run the filter function
    filter_fasta_dir(args.input_dir, args.output_fasta, args.patterns)


if __name__ == "__main__":
    main()
