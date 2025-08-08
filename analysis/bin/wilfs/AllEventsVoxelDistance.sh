#!/bin/bash

# Check arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_folder>"
    echo "Example: $0 /path/to/root_files/"
    echo ""
    echo "This script will:"
    echo "  - Find all .root files in the input folder"
    echo "  - Run DistanceAnalysis_3D on each file"
    echo "  - Generate event_XXXXXX_distanceGood.png files"
    exit 1
fi

INPUT_FOLDER="$1"

# Check if input folder exists
if [ ! -d "$INPUT_FOLDER" ]; then
    echo "Error: Input folder $INPUT_FOLDER not found"
    exit 1
fi

echo "Processing all ROOT files from folder: $INPUT_FOLDER"
echo ""

# Find all .root files in the input folder
ROOT_FILES=$(find "$INPUT_FOLDER" -name "*.root" -type f | sort)

if [ -z "$ROOT_FILES" ]; then
    echo "No .root files found in $INPUT_FOLDER"
    exit 1
fi

# Count total files
TOTAL_FILES=$(echo "$ROOT_FILES" | wc -l)
echo "Found $TOTAL_FILES ROOT files to process:"
echo "$ROOT_FILES"
echo ""

# Initialize counters
PROCESSED=0
SUCCESSFUL=0
FAILED=0

# Process each ROOT file
while IFS= read -r ROOT_FILE; do
    PROCESSED=$((PROCESSED + 1))
    BASENAME=$(basename "$ROOT_FILE")
    
    echo "[$PROCESSED/$TOTAL_FILES] Processing: $BASENAME"
    echo "  File: $ROOT_FILE"
    
    # Run the distance analysis
    ./VoxelDistanceFilter "$ROOT_FILE"
    
    # Check if the command succeeded
    if [ $? -eq 0 ]; then
        echo "  ✓ Successfully processed $BASENAME"
        SUCCESSFUL=$((SUCCESSFUL + 1))
    else
        echo "  ✗ Error processing $BASENAME"
        FAILED=$((FAILED + 1))
    fi
    
    echo "---"
done <<< "$ROOT_FILES"

echo ""
echo "Batch processing complete!"
echo "Summary:"
echo "  Total files: $TOTAL_FILES"
echo "  Successful: $SUCCESSFUL"
echo "  Failed: $FAILED"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "Warning: $FAILED files failed to process"
    exit 1
else
    echo ""
    echo "All files processed successfully!"
    echo "PNG files have been saved to their respective FitDistance_PNG directories"
fi