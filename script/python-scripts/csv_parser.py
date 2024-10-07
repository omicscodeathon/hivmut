import os
import pandas as pd
import argparse


def convert_txt_to_csv(directory, output_file):
    # Create an empty list to store DataFrames
    all_data = []

    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):  # Only process .txt files
            file_path = os.path.join(directory, filename)

            # Read the tab-delimited .txt file into a DataFrame
            df = pd.read_csv(file_path, sep="\t", index_col=False)

            # Append the DataFrame to the list
            all_data.append(df)
            print(f"Processed file: {filename}")

    # Concatenate all DataFrames into one
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)

        # Write the combined data to a CSV file
        combined_data.to_csv(output_file, index=False)
        print(f"Data successfully saved to {output_file}")
    else:
        print("No .txt files found in the specified directory.")


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Convert tab-delimited .txt files to a single .csv file."
    )

    # Add arguments
    parser.add_argument(
        "directory", type=str, help="Directory containing the .txt files."
    )
    parser.add_argument("output_file", type=str, help="Output CSV file name.")

    # Parse arguments
    args = parser.parse_args()

    # Run the conversion with the provided arguments
    convert_txt_to_csv(args.directory, args.output_file)
