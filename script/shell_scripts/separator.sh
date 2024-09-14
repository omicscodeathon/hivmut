#!/bin/bash

# Usage function to show help information
usage() {
    echo "Usage: $0 -d <directories> -o <output_dir> -p <pattern1> [-p <pattern2> ...] -n <batch_size>"
    echo "  -d  Directories to search (comma-separated if multiple)"
    echo "  -o  Output directory for the concatenated FASTA files"
    echo "  -p  Patterns to search for in the files (multiple allowed)"
    echo "  -n  Number of files per batch"
    exit 1
}

# Initialize variables
DIRECTORIES=()
PATTERNS=()
OUTPUT_DIR=""
BATCH_SIZE=0

# Parse command-line arguments
while getopts ":d:o:p:n:" opt; do
  case ${opt} in
    d )
      IFS=',' read -r -a DIRECTORIES <<< "$OPTARG"
      ;;
    o )
      OUTPUT_DIR="$OPTARG"
      ;;
    p )
      PATTERNS+=("$OPTARG")
      ;;
    n )
      BATCH_SIZE="$OPTARG"
      ;;
    \? )
      usage
      ;;
  esac
done

# Check if required arguments are provided
if [ -z "$DIRECTORIES" ] || [ -z "$OUTPUT_DIR" ] || [ "${#PATTERNS[@]}" -eq 0 ] || [ "$BATCH_SIZE" -le 0 ]; then
    usage
fi

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Temporary file to hold matched sequences before batching
TEMP_FILE="temp_filtered_sequences.fasta"
> "$TEMP_FILE"  # Clear temp file

# Counter to keep track of batch numbers and files processed
BATCH_COUNT=1
SEQUENCE_COUNT=0
TOTAL_FILES=0
MATCHED_FILES=0

# Function to process a directory
process_directory() {
  local DIR=$1
  local TOTAL_FILES_DIR=0
  for FILE in "$DIR"/*.fasta; do
    TOTAL_FILES=$((TOTAL_FILES + 1))
    TOTAL_FILES_DIR=$((TOTAL_FILES_DIR + 1))

    # Search for patterns
    MATCHED=true
    for PATTERN in "${PATTERNS[@]}"; do
      if ! grep -iq "$PATTERN" "$FILE"; then
        MATCHED=false
        break
      fi
    done

    if [ "$MATCHED" = true ]; then
      cat "$FILE" >> "$TEMP_FILE"
      SEQUENCE_COUNT=$((SEQUENCE_COUNT + 1))
      MATCHED_FILES=$((MATCHED_FILES + 1))

      # Display progress
      echo -ne "Processing: $TOTAL_FILES files checked, $MATCHED_FILES matched\r"

      # If we have reached the batch size, create a new batch file
      if [ "$SEQUENCE_COUNT" -eq "$BATCH_SIZE" ]; then
        mv "$TEMP_FILE" "$OUTPUT_DIR/batch_${BATCH_COUNT}.fasta"
        BATCH_COUNT=$((BATCH_COUNT + 1))
        SEQUENCE_COUNT=0
        > "$TEMP_FILE"  # Clear temp file for the next batch
      fi
    fi
  done
  echo -ne "\nFinished processing $TOTAL_FILES_DIR files in $DIR\n"
}

# Process all directories
for DIR in "${DIRECTORIES[@]}"; do
  process_directory "$DIR"
done

# If there are any remaining sequences that didn't make it to a full batch, save them
if [ "$SEQUENCE_COUNT" -gt 0 ]; then
  mv "$TEMP_FILE" "$OUTPUT_DIR/batch_${BATCH_COUNT}.fasta"
fi

# Clean up temporary file if it still exists
[ -f "$TEMP_FILE" ] && rm "$TEMP_FILE"

echo -e "\nProcessing complete. $MATCHED_FILES FASTA files were matched and saved in $OUTPUT_DIR."
