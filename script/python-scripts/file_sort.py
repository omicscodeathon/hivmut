"""
file_sort.py

This script sorts and moves FASTA files based on a list of accession IDs provided in a CSV file. 
It searches for these files in specified source directories and copies them to a designated destination directory.

Author: Halleluyah Darasimi Oludele
Date: September 2024

Usage:
    python file_sort.py --directories <source_dir1> <source_dir2> ... --destination <dest_dir> --csv-file <csv_file_path>

Arguments:
    --directories       List of source directories to search for FASTA files.
    --destination       Path to the destination directory where files will be copied.
    --csv-file          Path to the CSV file containing accession IDs.

Example:
    python file_sort.py --directories "../../data/hiv_fasta_files/" "../../data/Full Genome/" --destination "../../data/analysis-genome/" --csv-file "../../data/processed/hiv_lanl_analysis.csv"
"""

import os
import shutil
import pandas as pd
import argparse
import logging


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )


def read_accessions(csv_file_path):
    """Read accession IDs from a CSV file."""
    try:
        df = pd.read_csv(csv_file_path)
        return df["Accession"].tolist()
    except Exception as e:
        logging.error(f"Error reading CSV file {csv_file_path}: {e}")
        raise


def create_destination_directory(destination_path):
    """Create destination directory if it does not exist."""
    if not os.path.exists(destination_path):
        os.makedirs(destination_path)
        logging.info(f"Created destination directory {destination_path}")


def copy_files(file_list, directories, destination_path):
    """Copy files from source directories to the destination directory."""
    for file in file_list:
        file_copied = False
        for directory in directories:
            file_path = os.path.join(directory, file)
            if os.path.isfile(file_path):
                shutil.copy(file_path, destination_path)
                logging.info(f"Copied {file} from {directory} to {destination_path}")
                file_copied = True
                break
        if not file_copied:
            logging.warning(f"File {file} not found in any of the directories")


def main(args):
    """Main function to execute the file sorting program."""
    setup_logging()

    # Read accession list from CSV
    accession_list = read_accessions(args.csv_file)

    # Prepare list of files to copy
    file_list = [accession + ".fasta" for accession in accession_list]

    # Create destination directory
    create_destination_directory(args.destination)

    # Copy files from source directories to destination
    copy_files(file_list, args.directories, args.destination)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sort and move FASTA files based on accession list."
    )
    parser.add_argument(
        "--directories",
        nargs="+",
        required=True,
        help="List of source directories to search for FASTA files.",
    )
    parser.add_argument(
        "--destination",
        required=True,
        help="Destination directory where files will be copied.",
    )
    parser.add_argument(
        "--csv-file",
        required=True,
        help="Path to the CSV file containing accession IDs.",
    )

    args = parser.parse_args()
    main(args)
