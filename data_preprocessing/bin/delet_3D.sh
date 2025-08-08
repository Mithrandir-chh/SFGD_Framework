#!/bin/bash

# Script to delete all folders containing "3D_Display_root" in their names
# Usage: ./delete_3d_display_folders.sh /path/to/main/folder

# Check if path argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No path provided"
    echo "Usage: $0 /path/to/main/folder"
    exit 1
fi

# Check if the provided path exists
if [ ! -d "$1" ]; then
    echo "Error: Directory '$1' does not exist"
    exit 1
fi

SEARCH_PATH="$1"
SEARCH_PATTERN="*Slot*"

echo "Searching for folders containing '3D_Display_root' in: $SEARCH_PATH"
echo

# Find and display folders that will be deleted
# folders_to_delete=$(find "$SEARCH_PATH" -type d -name "$SEARCH_PATTERN" 2>/dev/null)

# If deleting files
folders_to_delete=$(find "$SEARCH_PATH" -type f -name "$SEARCH_PATTERN" 2>/dev/null)

if [ -z "$folders_to_delete" ]; then
    echo "No folders containing '3D_Display_root' found in $SEARCH_PATH"
    exit 0
fi

echo "Found the following folders to delete:"
echo "$folders_to_delete"
echo

read -p "Are you sure you want to delete these folders and all their contents? (y/N): " confirm

case $confirm in
    [yY]|[yY][eE][sS])
        echo "Deleting folders..."
        
        # Delete each folder
        while IFS= read -r folder; do
            # if [ -d "$folder" ]; then
            if [ -f "$folder" ]; then   # if deleting files         
                echo "Deleting: $folder"
                # rm -rf "$folder"
                rm "$folder" # if deleting files
                if [ $? -eq 0 ]; then
                    echo "✓ Successfully deleted: $folder"
                else
                    echo "✗ Failed to delete: $folder"
                fi
            fi
        done <<< "$folders_to_delete"
        
        echo "Deletion process completed."
        ;;
    *)
        echo "Operation cancelled."
        exit 0
        ;;
esac