#!/bin/bash

# --- Configuration ---
# Set the maximum number of subfolders to process in parallel.
MAX_JOBS=20
# Set the final destination directory for the SingleTrack results.
Beta_DESTINATION_DIR="/media/disk_b/standard_software/sfgd_framework/cosmic_data_2025/che_beta1/TotalBetaNarrow"
Gamma_DESTINATION_DIR="/media/disk_b/standard_software/sfgd_framework/cosmic_data_2025/che_beta1/Ditch"

# --- Function Definitions ---

# This function contains all the logic to process a SINGLE subfolder.
process_subfolder() {
    local subfolder=$1
    local beta_dest_dir=$2 # Receive the destination directory as an argument
    local gamma_dest_dir=$3
    local subfolder_name=$(basename "$subfolder")

    echo "[$(date '+%H:%M:%S')] Starting: $subfolder_name"

    local event_files_found=false
    # Loop through event files in the subfolder
    for f in "${subfolder}/"*events.root; do
        if [ ! -e "$f" ]; then continue; fi # Skip if glob finds no files

        event_files_found=true
        echo "[$(date '+%H:%M:%S')]   Displaying: $(basename "$f") in $subfolder_name"
        if ./3D_EventGamma "$f"; then
            echo "[$(date '+%H:%M:%S')]   ✓ MuonDecay successful for $(basename "$f")"
        else
            echo "[$(date '+%H:%M:%S')]   ✗ MuonDecay FAILED for $(basename "$f")"
            continue # Skip to the next event file if this one failed
        fi
    done

    # --- Copy MuonDecay contents after all processing is complete for this subfolder ---
    local MuonDecay_copied=false
    for d in "${subfolder}/"*BetaEvents*; do
        if [ -d "$d" ]; then
            # --- NEW: Check if the directory has any files before trying to copy ---
            if [ -z "$(ls -A "$d")" ]; then
                echo "[$(date '+%H:%M:%S')]   - Skipping copy from $(basename "$d"), directory is empty."
                continue
            fi
            
            echo "[$(date '+%H:%M:%S')]   Copying results from $(basename "$d") to destination"
            # Use -n to not overwrite existing files, which is safer in parallel
            # REMOVED 2>/dev/null to make errors visible
            if cp -r -n "$d"/. "$beta_dest_dir"/; then
                echo "[$(date '+%H:%M:%S')]   ✓ Copy successful from $(basename "$d")"
                MuonDecay_copied=true
            else
                # This now indicates a real failure (e.g., permissions)
                echo "[$(date '+%H:%M:%S')]   ✗ Copy FAILED from $(basename "$d") (exit code: $?)"
            fi
        fi
    done

    for d in "${subfolder}/"*GammaEvents*; do
        if [ -d "$d" ]; then
            # --- NEW: Check if the directory has any files before trying to copy ---
            if [ -z "$(ls -A "$d")" ]; then
                echo "[$(date '+%H:%M:%S')]   - Skipping copy from $(basename "$d"), directory is empty."
                continue
            fi
            
            echo "[$(date '+%H:%M:%S')]   Copying results from $(basename "$d") to destination"
            # Use -n to not overwrite existing files, which is safer in parallel
            # REMOVED 2>/dev/null to make errors visible
            if cp -r -n "$d"/. "$gamma_dest_dir"/; then
                echo "[$(date '+%H:%M:%S')]   ✓ Copy successful from $(basename "$d")"
                MuonDecay_copied=true
            else
                # This now indicates a real failure (e.g., permissions)
                echo "[$(date '+%H:%M:%S')]   ✗ Copy FAILED from $(basename "$d") (exit code: $?)"
            fi
        fi
    done

    if [ "$MuonDecay_copied" = false ]; then
        echo "[$(date '+%H:%M:%S')]   → No '*MuonDecay' directory with content found in $subfolder_name to copy"
    fi

    echo "[$(date '+%H:%M:%S')] ✓ Finished: $subfolder_name"
}

# Export the function so it's available to the parallel sub-shells
export -f process_subfolder

# --- Main Script Logic ---

# Check if at least one folder path is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 /path/to/folder1/ [/path/to/folder2/] ..."
    echo "Example: $0 /media/disk_b/.../che_side0/"
    exit 1
fi

# Ensure the destination directory exists
mkdir -p "$Beta_DESTINATION_DIR"
echo "Starting processing with up to $MAX_JOBS parallel jobs."
echo "Results will be copied to: $Beta_DESTINATION_DIR"
echo "========================================================"

# First, collect all valid subfolders from all input directories
all_subfolders=()
for FOLDER_PATH in "$@"; do
    FOLDER_PATH="${FOLDER_PATH%/}"
    if [ ! -d "$FOLDER_PATH" ]; then
        echo "Warning: Input folder '$FOLDER_PATH' does not exist. Skipping..."
        continue
    fi
    while IFS= read -r -d '' sub; do
        all_subfolders+=("$sub")
    done < <(find "$FOLDER_PATH" -mindepth 1 -maxdepth 1 -type d -print0)
done

echo "Found ${#all_subfolders[@]} subfolders to process."

# Now, loop through the collected subfolders and process them in parallel
job_count=0
for subfolder in "${all_subfolders[@]}"; do
    # Corrected job control logic
    while (( $(jobs -p | wc -l) >= MAX_JOBS )); do
        sleep 0.2
    done
    
    # Run the processing function in the background, passing the destination directory
    process_subfolder "$subfolder" "$Beta_DESTINATION_DIR" "$Gamma_DESTINATION_DIR"&
    job_count=$((job_count + 1))
done

# Wait for all background jobs to complete
echo "========================================================"
echo "All $job_count jobs launched. Waiting for completion..."
wait
echo "All folders processed."
