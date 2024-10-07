import argparse
import csv


def split_csv(input_file, output_prefix, chunk_size=70000):
    """Splits a CSV file into smaller files with a specified chunk size.

    Args:
        input_file (str): The path to the input CSV file.
        output_prefix (str): The prefix for the output file names.
        chunk_size (int): The maximum number of rows per output file.
    """

    with open(input_file, "r") as f:
        reader = csv.reader(f)
        header = next(reader)  # Read the header row
        rows = list(reader)

        num_chunks = (len(rows) - 1) // chunk_size + 1
        for i in range(num_chunks):
            start_index = i * chunk_size
            end_index = min((i + 1) * chunk_size, len(rows))
            chunk = rows[start_index:end_index]

            output_file = f"{output_prefix}_{i + 1}.csv"
            with open(output_file, "w", newline="") as output:
                writer = csv.writer(output)
                writer.writerow(header)
                writer.writerows(chunk)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a CSV file into smaller chunks."
    )
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_prefix", help="Prefix for the output file names")
    parser.add_argument(
        "-c",
        "--chunk_size",
        type=int,
        default=70000,
        help="Maximum number of rows per output file",
    )

    args = parser.parse_args()

    split_csv(args.input_file, args.output_prefix, args.chunk_size)



# def extract_and_convert(message):
#     # Regex to detect missing genes (e.g., PR, RT, IN)
#     missing_genes = list(set(re.findall(r"\bPR\b|\bRT\b|\bIN\b", message)))

#     # Regex to detect drug resistance positions
#     resistance_positions = list(set(re.findall(r"RT\s\d+", message)))

#     # Look for specific issues (frameshift, missing sequences, etc.)
#     issues = []
#     if "frameshift" in message.lower():
#         issues.append("frameshift")
#     if "refuse to process" in message.lower():
#         issues.append("refuse to process")
#     if "not sequenced or aligned" in message.lower():
#         issues.append("missing sequence")

#     # Convert lists to comma-separated strings or NaN for empty lists
#     missing_genes_str = ", ".join(missing_genes) if missing_genes else np.nan
#     resistance_positions_str = (
#         ", ".join(resistance_positions) if resistance_positions else np.nan
#     )
#     issues_str = ", ".join(issues) if issues else np.nan

#     return pd.Series([missing_genes_str, resistance_positions_str, issues_str])


# # Apply the combined function and expand the results into new columns
# validation_df[["Missing Genes", "Drug Resistance Positions", "Issues"]] = validation_df[
#     "Message"
# ].apply(extract_and_convert)

# # Display the extracted info
# validation_df.head()
