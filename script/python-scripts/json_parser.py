import json
import csv

# Load your JSON file
with open("your_json_file.json", "r") as file:
    data = json.load(file)

# Define a CSV file and the headers (columns) based on the schema fields you want to extract
csv_file = "output.csv"
csv_columns = [
    "inputSequence.header",
    "inputSequence.SHA512",
    "strain.name",
    "subtypeText",
    "alignedGeneSequences.gene.name",
    "alignedGeneSequences.gene.length",
    "alignedGeneSequences.firstAA",
    "alignedGeneSequences.lastAA",
    "mutation.position",
    "mutation.AAs",
    "mutation.isInsertion",
    "mutation.isDeletion",
    "mutation.isApobecMutation",
    "drugResistance.gene.name",
    "drugResistance.drugScores.drug.name",
    "drugResistance.drugScores.score",
    "drugResistance.drugScores.level",
]

# Open CSV file for writing
with open(csv_file, mode="w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()

    # Iterate over each data entry in the JSON
    for entry in data:
        row = {}

        # Extract input sequence data
        row["inputSequence.header"] = entry["inputSequence"]["header"]
        row["inputSequence.SHA512"] = entry["inputSequence"]["SHA512"]

        # Extract strain and subtype data
        row["strain.name"] = entry["strain"]["name"]
        row["subtypeText"] = entry["subtypeText"]

        # Iterate over aligned gene sequences
        for gene_seq in entry["alignedGeneSequences"]:
            row["alignedGeneSequences.gene.name"] = gene_seq["gene"]["name"]
            row["alignedGeneSequences.gene.length"] = gene_seq["gene"]["length"]
            row["alignedGeneSequences.firstAA"] = gene_seq["firstAA"]
            row["alignedGeneSequences.lastAA"] = gene_seq["lastAA"]

            # Extract mutations for each gene sequence
            for mutation in gene_seq["mutations"]:
                row["mutation.position"] = mutation["position"]
                row["mutation.AAs"] = mutation["AAs"]
                row["mutation.isInsertion"] = mutation["isInsertion"]
                row["mutation.isDeletion"] = mutation["isDeletion"]
                row["mutation.isApobecMutation"] = mutation["isApobecMutation"]

                # Write the row after extracting all mutation fields
                writer.writerow(row)

        # Extract drug resistance data
        for drug_resistance in entry["drugResistance"]:
            row["drugResistance.gene.name"] = drug_resistance["gene"]["name"]
            for drug_score in drug_resistance["drugScores"]:
                row["drugResistance.drugScores.drug.name"] = drug_score["drug"]["name"]
                row["drugResistance.drugScores.score"] = drug_score["score"]
                row["drugResistance.drugScores.level"] = drug_score["level"]

                # Write the row with drug resistance data
                writer.writerow(row)
