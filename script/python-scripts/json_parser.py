import json
import csv
import os
import argparse
from tqdm import tqdm


def extract_general_mutations(entry):
    """Extracts general mutation information from the SierraPy JSON entry."""
    mutation_data = []

    for aligned_gene in entry.get("alignedGeneSequences", []):
        gene_name = aligned_gene.get("gene", {}).get("name", None)
        mutations = aligned_gene.get("mutations", [])

        for mutation in mutations:
            mutation_data.append(
                {
                    "Input Sequence": entry.get("inputSequence", {}).get(
                        "header", None
                    ),
                    "SHA512": entry.get("inputSequence", {}).get("SHA512", None),
                    "Strain Name": entry.get("strain", {}).get("name", None),
                    "Subtype": entry.get("subtypeText", None),
                    "Gene Name": gene_name,
                    "Mutation Position": mutation.get("position", None),
                    "Amino Acids": mutation.get("AAs", None),
                    "Primary Type": mutation.get("primaryType", None),
                    "Consensus": mutation.get("consensus", None),
                    "Is Insertion": mutation.get("isInsertion", None),
                    "Is Deletion": mutation.get("isDeletion", None),
                    "Is Indel": mutation.get("isIndel", None),
                    "Is Apobec Mutation": mutation.get("isApobecMutation", None),
                    "Is Apobec DRM": mutation.get("isApobecDRM", None),
                    "Has Stop Codon": mutation.get("hasStop", None),
                    "Is Unusual": mutation.get("isUnusual", None),
                    "Is SDRM": mutation.get("isSDRM", None),
                    "Mutation Text": mutation.get("text", None),
                }
            )
    return mutation_data


def extract_drm_info(entry):
    """Extracts drug resistance mutation (DRM) information from the SierraPy JSON entry."""
    drm_data = []

    for drug_resistance in entry.get("drugResistance", []):
        version_text = drug_resistance.get("version", {}).get("text", None)
        publish_date = drug_resistance.get("version", {}).get("publishDate", None)
        gene_name = drug_resistance.get("gene", {}).get("name", None)

        for drug_score in drug_resistance.get("drugScores", []):
            for partial_score in drug_score.get("partialScores", []):
                for mutation in partial_score.get("mutations", []):
                    drm_data.append(
                        {
                            "Input Sequence": entry.get("inputSequence", {}).get(
                                "header", None
                            ),
                            "Gene Name": gene_name,
                            "Drug Class": drug_score.get("drugClass", {}).get(
                                "name", None
                            ),
                            "Drug": drug_score.get("drug", {}).get("displayAbbr", None),
                            "Drug Full Name": drug_score.get("drug", {}).get(
                                "name", None
                            ),
                            "Score": drug_score.get("score", None),
                            "Level": drug_score.get("level", None),
                            "Version": version_text,
                            "Publish Date": publish_date,
                            "Mutation Text": mutation.get("text", None),
                            "Mutation Primary Type": mutation.get("primaryType", None),
                            "Comments": mutation.get("comments", None),
                        }
                    )
    return drm_data


def extract_gene_info(entry):
    """Extracts gene information from the SierraPy JSON entry."""
    gene_data = []

    for aligned_gene in entry.get("alignedGeneSequences", []):
        gene_data.append(
            {
                "Input Sequence": entry.get("inputSequence", {}).get("header", None),
                "Strain Name": entry.get("strain", {}).get("name", None),
                "Subtype": entry.get("subtypeText", None),
                "Gene Name": aligned_gene.get("gene", {}).get("name", None),
                "First AA": aligned_gene.get("firstAA", None),
                "Last AA": aligned_gene.get("lastAA", None),
                "First NA": aligned_gene.get("firstNA", None),
                "Last NA": aligned_gene.get("lastNA", None),
                "Gene Length": aligned_gene.get("gene", {}).get("length", None),
                "Match Percentage": aligned_gene.get("matchPcnt", None),
            }
        )
    return gene_data


def extract_pretty_pairwise(entry):
    """Extracts the pretty pairwise alignment information."""
    pairwise_data = []

    for aligned_gene in entry.get("alignedGeneSequences", []):
        pairwise = aligned_gene.get("prettyPairwise", {})
        pairwise_data.append(
            {
                "Input Sequence": entry.get("inputSequence", {}).get("header", None),
                "Gene Name": aligned_gene.get("gene", {}).get("name", None),
                "Position Line": "; ".join(pairwise.get("positionLine", [])) or None,
                "Reference AA Line": "; ".join(pairwise.get("refAALine", [])) or None,
                "Aligned NA Line": "; ".join(pairwise.get("alignedNAsLine", []))
                or None,
                "Mutation Line": "; ".join(pairwise.get("mutationLine", [])) or None,
            }
        )
    return pairwise_data


def extract_validation_results(entry):
    """Extracts validation results from the SierraPy JSON entry."""
    validation_data = []

    for result in entry.get("validationResults", []):
        validation_data.append(
            {
                "Input Sequence": entry.get("inputSequence", {}).get("header", None),
                "Level": result.get("level", None),
                "Message": result.get("message", None),
            }
        )
    return validation_data


def parse_sierra_json(json_data, accessions_set):
    """Parses SierraPy JSON data and extracts relevant information, accumulating accessions."""
    general_mutation_data = []
    drm_info_data = []
    gene_info_data = []
    pairwise_data = []
    validation_data = []

    for entry in json_data:
        # Add the input sequence header to the accessions set
        accessions_set.add(entry["inputSequence"]["header"])

        general_mutation_data.extend(extract_general_mutations(entry))
        drm_info_data.extend(extract_drm_info(entry))
        gene_info_data.extend(extract_gene_info(entry))
        pairwise_data.extend(extract_pretty_pairwise(entry))
        validation_data.extend(extract_validation_results(entry))

    return (
        general_mutation_data,
        drm_info_data,
        gene_info_data,
        pairwise_data,
        validation_data,
    )


def save_to_csv(parsed_data, output_file, headers):
    """Saves the parsed data to a CSV file."""
    with open(output_file, mode="a", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=headers)

        # Write header only if file is empty
        if file.tell() == 0:
            writer.writeheader()

        writer.writerows(parsed_data)


def save_accessions(accessions_set, output_file):
    """Saves all unique accessions to a text file."""
    with open(output_file, "w") as file:
        for accession in sorted(accessions_set):
            file.write(f"{accession}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Parse multiple SierraPy JSON files and save combined data to CSV."
    )
    parser.add_argument(
        "input_dir", help="Path to the input directory containing SierraPy JSON files."
    )
    parser.add_argument("output_dir", help="Directory to save the output CSV files.")

    args = parser.parse_args()

    # Prepare output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Set of all accessions found
    accessions_set = set()

    # Initialize combined data
    combined_general_mutation_data = []
    combined_drm_info_data = []
    combined_gene_info_data = []
    combined_pairwise_data = []
    combined_validation_data = []

    # List all JSON files in the input directory
    json_files = [f for f in os.listdir(args.input_dir) if f.endswith(".json")]

    # Process each JSON file with tqdm progress bar
    with tqdm(total=len(json_files), desc="Processing JSON files", unit="file") as pbar:
        for json_file in json_files:
            with open(os.path.join(args.input_dir, json_file), "r") as f:
                json_data = json.load(f)
                (
                    general_mutation_data,
                    drm_info_data,
                    gene_info_data,
                    pairwise_data,
                    validation_data,
                ) = parse_sierra_json(json_data, accessions_set)

                # Append the parsed data to the combined lists
                combined_general_mutation_data.extend(general_mutation_data)
                combined_drm_info_data.extend(drm_info_data)
                combined_gene_info_data.extend(gene_info_data)
                combined_pairwise_data.extend(pairwise_data)
                combined_validation_data.extend(validation_data)

            # Update progress bar after processing each file
            pbar.update(1)

    # Save combined data to CSV files
    save_to_csv(
        combined_general_mutation_data,
        os.path.join(args.output_dir, "general_mutation_data.csv"),
        combined_general_mutation_data[0].keys()
        if combined_general_mutation_data
        else [],
    )
    save_to_csv(
        combined_drm_info_data,
        os.path.join(args.output_dir, "drug_resistance_info.csv"),
        combined_drm_info_data[0].keys() if combined_drm_info_data else [],
    )
    save_to_csv(
        combined_gene_info_data,
        os.path.join(args.output_dir, "gene_info_data.csv"),
        combined_gene_info_data[0].keys() if combined_gene_info_data else [],
    )
    save_to_csv(
        combined_pairwise_data,
        os.path.join(args.output_dir, "pretty_pairwise_data.csv"),
        combined_pairwise_data[0].keys() if combined_pairwise_data else [],
    )
    save_to_csv(
        combined_validation_data,
        os.path.join(args.output_dir, "validation_data.csv"),
        combined_validation_data[0].keys() if combined_validation_data else [],
    )

    # Save all unique accessions to a text file
    save_accessions(
        accessions_set, os.path.join(args.output_dir, "unique_accessions.txt")
    )


if __name__ == "__main__":
    main()
