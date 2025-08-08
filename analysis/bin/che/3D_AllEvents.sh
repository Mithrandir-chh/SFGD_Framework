#!/bin/bash

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file.root> <png_folder_path>"
    echo "Example: $0 /path/to/events.root /path/to/event_pngs/"
    exit 1
fi

ROOT_FILE="$1"
PNG_FOLDER="$2"

# Check if ROOT file exists
if [ ! -f "$ROOT_FILE" ]; then
    echo "Error: ROOT file $ROOT_FILE not found"
    exit 1
fi

# Check if PNG folder exists
if [ ! -d "$PNG_FOLDER" ]; then
    echo "Error: PNG folder $PNG_FOLDER not found"
    exit 1
fi

echo "Processing all events from PNG folder: $PNG_FOLDER"
echo "Using ROOT file: $ROOT_FILE"
echo ""

# Extract event numbers from PNG filenames and sort them (FIXED)
EVENT_NUMBERS=$(ls "$PNG_FOLDER"/event_*.png 2>/dev/null | \
    sed 's/.*event_\([0-9]\+\)\(_.*\)\?\.png/\1/' | \
    sort -n)

if [ -z "$EVENT_NUMBERS" ]; then
    echo "No event PNG files found in $PNG_FOLDER"
    exit 1
fi

echo "Found events: $EVENT_NUMBERS"
echo ""

# Process each event
for EVENT_NUM in $EVENT_NUMBERS; do
    echo "Processing event $EVENT_NUM..."
    ./3D_EventTrajectory "$ROOT_FILE" "$EVENT_NUM"
    
    # Check if the command succeeded
    if [ $? -ne 0 ]; then
        echo "Error processing event $EVENT_NUM"
        exit 1
    fi
    echo "---"
done

echo "All events processed successfully!"