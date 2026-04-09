#!/bin/bash

# Script: extract_fasta_region.sh
# Usage: ./extract_fasta_region.sh <fasta_file> <sequence_id> <start> <end>
# Example: ./extract_fasta_region.sh genome.fasta xyz 1000 2000

# Check if correct number of arguments provided
if [ $# -ne 4 ]; then
    echo "Usage: $0 <fasta_file> <sequence_id> <start> <end>"
    echo "Example: $0 genome.fasta chr1 1000 2000"
    echo ""
    echo "Arguments:"
    echo "  fasta_file    : Input FASTA file"
    echo "  sequence_id   : Sequence identifier (without '>')"
    echo "  start         : Start coordinate (1-based, inclusive)"
    echo "  end           : End coordinate (1-based, inclusive)"
    exit 1
fi

FASTA_FILE="$1"
SEQ_ID="$2"
START="$3"
END="$4"

# Validate that start and end are positive integers
if ! [[ "$START" =~ ^[0-9]+$ ]] || ! [[ "$END" =~ ^[0-9]+$ ]]; then
    echo "Error: Start and end must be positive integers"
    exit 1
fi

# Check if start is less than or equal to end
if [ "$START" -gt "$END" ]; then
    echo "Error: Start coordinate ($START) must be less than or equal to end coordinate ($END)"
    exit 1
fi

# Check if FASTA file exists
if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: FASTA file '$FASTA_FILE' not found"
    exit 1
fi

# Calculate sequence length
SEQ_LENGTH=$(awk -v id="$SEQ_ID" '
    /^>/ {
        if (match($0, id)) found=1; else found=0
        next
    }
    found {seq = seq $0}
    END {print length(seq)}
' "$FASTA_FILE")

# Check if sequence ID was found
if [ -z "$SEQ_LENGTH" ] || [ "$SEQ_LENGTH" -eq 0 ]; then
    echo "Error: Sequence ID '$SEQ_ID' not found in FASTA file"
    exit 1
fi

# Check if coordinates are within bounds
if [ "$START" -gt "$SEQ_LENGTH" ] || [ "$END" -gt "$SEQ_LENGTH" ]; then
    echo "Error: Coordinates exceed sequence length ($SEQ_LENGTH)"
    echo "  Requested: $START-$END"
    echo "  Available: 1-$SEQ_LENGTH"
    exit 1
fi

# Extract the region
echo "Extracting region $SEQ_ID:$START-$END from $FASTA_FILE..."
echo ">${SEQ_ID}_${START}_${END}"

# Method 1: Using awk (more memory efficient for large sequences)
awk -v id="$SEQ_ID" -v start="$START" -v end="$END" '
    /^>/ {
        if (match($0, id)) found=1; else found=0
        if (found) seq=""
        next
    }
    found {seq = seq $0}
    END {
        if (found) {
            region = substr(seq, start, end-start+1)
            # Format output in lines of 60 characters
            for(i=1; i<=length(region); i+=60) {
                print substr(region, i, 60)
            }
        }
    }
' "$FASTA_FILE"

# Alternative Method 2: Using bioawk (if available, handles line breaks automatically)
# Uncomment the following lines if you have bioawk installed
# if command -v bioawk &> /dev/null; then
#     echo ">${SEQ_ID}_${START}_${END}"
#     bioawk -c fastx -v seq="$SEQ_ID" -v s="$START" -v e="$END" \
#         '$name == seq {print substr($seq, s, e-s+1)}' "$FASTA_FILE" | \
#         fold -w 60
# fi

echo ""
echo "Extraction complete!"
