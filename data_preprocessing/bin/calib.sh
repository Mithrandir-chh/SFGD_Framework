#!/bin/bash

# --- Configuration ---
# Set the maximum number of subfolders to process in parallel.
# A good starting point is the number of CPU cores you have, or slightly less.
MAX_JOBS=20

# --- Function Definitions ---

# This function contains all the logic to process a SINGLE subfolder.
# It will be run in the background for each subfolder.
process_subfolder() {
    local subfolder=$1
    local subfolder_name=$(basename "$subfolder")

    echo "[$(date '+%H:%M:%S')] Starting: $subfolder_name"

    # --- Step 1: Run Calibration if needed ---
    # Check if any file ending with *calib.root exists to decide whether to skip.
    local calib_exists=false
    for existing_calib in "${subfolder}/"*calib.root; do
        if [ -e "$existing_calib" ]; then
            calib_exists=true
            break # Found one, no need to check further
        fi
    done

    if [ "$calib_exists" = true ]; then
        echo "[$(date '+%H:%M:%S')] ⏭ Skipping Calibration for $subfolder_name (calib.root file already exists)"
    else
        local raw_files_found=false
        for f in "${subfolder}/"*raw.root; do
            if [ -e "$f" ]; then
                raw_files_found=true
                echo "[$(date '+%H:%M:%S')]   Calibrating: $(basename "$f") in $subfolder_name"
                # Run the Calibration command
                if ./Calibration "$f" "55" "50"; then
                    echo "[$(date '+%H:%M:%S')]   ✓ Calib successful for $(basename "$f")"
                else
                    echo "[$(date '+%H:%M:%S')]   ✗ Calib FAILED for $(basename "$f") (exit code: $?)"
                fi
            fi
        done
        if [ "$raw_files_found" = false ]; then
            echo "[$(date '+%H:%M:%S')]   → No raw.root files found in $subfolder_name"
        fi
    fi

    # --- Step 2: Run EventStructure ---
    local calib_files_found=false
    for f in "${subfolder}/"*calib.root; do
        if [ -e "$f" ]; then
            calib_files_found=true
            echo "[$(date '+%H:%M:%S')]   Structuring: $(basename "$f") in $subfolder_name"
            # Run the EventStructure command
            if ./EventStructure_12free "$f"; then
                echo "[$(date '+%H:%M:%S')]   ✓ EventStructure successful for $(basename "$f")"
            else
                echo "[$(date '+%H:%M:%S')]   ✗ EventStructure FAILED for $(basename "$f") (exit code: $?)"
            fi
        fi
    done

    if [ "$calib_files_found" = false ]; then
        echo "[$(date '+%H:%M:%S')]   → No calib.root files found in $subfolder_name to structure"
    fi

    echo "[$(date '+%H:%M:%S')] ✓ Finished: $subfolder_name"
}

# Export the function so it's available to the parallel sub-shells
export -f process_subfolder

# --- Main Script Logic ---

# Check if at least one folder path is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 /path/to/folder1/ [/path/to/folder2/] ..."
    echo "Example: $0 /media/disk_b/.../che_side0/ /media/disk_b/.../che_side1/"
    exit 1
fi

echo "Starting processing with up to $MAX_JOBS parallel jobs..."
echo "========================================================"

# First, collect all valid subfolders from all input directories into a single list
all_subfolders=()
for FOLDER_PATH in "$@"; do
    FOLDER_PATH="${FOLDER_PATH%/}" # Remove trailing slash
    if [ ! -d "$FOLDER_PATH" ]; then
        echo "Warning: Input folder '$FOLDER_PATH' does not exist. Skipping..."
        continue
    fi
    # Use find to get all first-level subdirectories
    while IFS= read -r -d '' sub; do
        all_subfolders+=("$sub")
    done < <(find "$FOLDER_PATH" -mindepth 1 -maxdepth 1 -type d -print0)
done

echo "Found ${#all_subfolders[@]} subfolders to process."

# Now, loop through the collected subfolders and process them in parallel
job_count=0
for subfolder in "${all_subfolders[@]}"; do
    # Wait if we have reached the maximum number of parallel jobs
    while (( $(jobs -p | wc -l) >= MAX_JOBS )); do
        sleep 0.2
    done
    
    # Run the processing function for the subfolder in the background
    process_subfolder "$subfolder" &
    job_count=$((job_count + 1))
done

# Wait for all background jobs to complete before exiting the script
echo "========================================================"
echo "All $job_count jobs launched. Waiting for completion..."
wait
echo "All folders processed."
