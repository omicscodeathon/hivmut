import argparse
import csv


def convert_txt_to_csv(input_path, output_path):
    with open(input_path, "r") as txt_file:
        lines = txt_file.readlines()  # Skip the first two lines

    with open(output_path, "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        for line in lines:
            csv_writer.writerow(
                line.strip().split("\t")
            )  # Assuming tab-separated values in the txt file


def main():
    parser = argparse.ArgumentParser(
        description="Convert a txt file to a csv file by removing the first two lines."
    )
    parser.add_argument("input_path", type=str, help="The path to the input txt file.")
    parser.add_argument(
        "output_path", type=str, help="The path to the output csv file."
    )

    args = parser.parse_args()
    convert_txt_to_csv(args.input_path, args.output_path)


if __name__ == "__main__":
    main()
