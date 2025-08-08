#!/bin/bash

# --- Configuration ---
MAX_JOBS=20
DESTINATION_DIR="/media/disk_b/standard_software/sfgd_framework/cosmic_data_2025/che_side0/Total/SingleTrack_000_to_017"

# Track PIDs for better cleanup
declare -a PIDS=()

# Cleanup function
cleanup() {
    echo "Cleaning up background processes..."
    for pid in "${PIDS[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            echo "Killing PID: $pid"
            kill "$pid" 2>/dev/null
        fi
    done
    wait
    exit
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM EXIT

# --- Function Definitions ---

process_subfolder() {
    local subfolder=$1
    local dest_dir=$2
    local subfolder_name=$(basename "$subfolder")

    echo "[$(date '+%H:%M:%S')] Starting: $subfolder_name"

    local event_files_found=false
    for f in "${subfolder}/"*events.root; do
        if [ ! -e "$f" ]; then continue; fi
        event_files_found=true

        local filtered_dir_exists=false
        for d_check in "${subfolder}/"*_Filtered; do
            if [ -d "$d_check" ]; then
                filtered_dir_exists=true
                break
            fi
        done

        if [ "$filtered_dir_exists" = true ]; then
            echo "[$(date '+%H:%M:%S')]   ⏭ Skipping EventDisplayLite for $(basename "$f"), Filtered directory already exists."
        else
            echo "[$(date '+%H:%M:%S')]   Displaying: $(basename "$f") in $subfolder_name"
            # ADD TIMEOUT to prevent hanging
            if timeout 300 ./EventDisplayLite -i "$f"; then
                echo "[$(date '+%H:%M:%S')]   ✓ EventDisplayLite successful for $(basename "$f")"
            else
                local exit_code=$?
                if [ $exit_code -eq 124 ]; then
                    echo "[$(date '+%H:%M:%S')]   ✗ EventDisplayLite TIMED OUT for $(basename "$f")"
                else
                    echo "[$(date '+%H:%M:%S')]   ✗ EventDisplayLite FAILED for $(basename "$f") (exit code: $exit_code)"
                fi
                continue
            fi
        fi

        local filtered_dir_found_for_3d=false
        for d in "${subfolder}/"*_Filtered; do
            if [ -d "$d" ]; then
                filtered_dir_found_for_3d=true
                echo "[$(date '+%H:%M:%S')]   Running 3D Script for $(basename "$d")"
                # ADD TIMEOUT here too
                if timeout 600 ./3D_AllEventsSide.sh "$f" "$d"; then
                    echo "[$(date '+%H:%M:%S')]   ✓ 3D_AllEventsSide.sh successful"
                else
                    local exit_code=$?
                    if [ $exit_code -eq 124 ]; then
                        echo "[$(date '+%H:%M:%S')]   ✗ 3D_AllEventsSide.sh TIMED OUT"
                    else
                        echo "[$(date '+%H:%M:%S')]   ✗ 3D_AllEventsSide.sh FAILED (exit code: $exit_code)"
                    fi
                fi
            fi
        done
        if [ "$filtered_dir_found_for_3d" = false ]; then
            echo "[$(date '+%H:%M:%S')]   → No '*Filtered' directory found to run 3D script for $(basename "$f")"
        fi
    done

    if [ "$event_files_found" = false ]; then
        echo "[$(date '+%H:%M:%S')]   → No *events.root files found in $subfolder_name"
    fi

    local singletrack_copied=false
    for d in "${subfolder}/"*SingleTrack_strict; do
        if [ -d "$d" ]; then
            if [ -z "$(ls -A "$d")" ]; then
                echo "[$(date '+%H:%M:%S')]   - Skipping copy from $(basename "$d"), directory is empty."
                continue
            fi
            
            echo "[$(date '+%H:%M:%S')]   Copying results from $(basename "$d") to destination"
            if cp -r -n "$d"/. "$dest_dir"/; then
                echo "[$(date '+%H:%M:%S')]   ✓ Copy successful from $(basename "$d")"
                singletrack_copied=true
            else
                echo "[$(date '+%H:%M:%S')]   ✗ Copy FAILED from $(basename "$d") (exit code: $?)"
            fi
        fi
    done

    if [ "$singletrack_copied" = false ]; then
        echo "[$(date '+%H:%M:%S')]   → No '*SingleTrack_strict' directory with content found in $subfolder_name to copy"
    fi

    echo "[$(date '+%H:%M:%S')] ✓ Finished: $subfolder_name"
}

export -f process_subfolder

# --- Main Script Logic ---

# Optional: Set total script timeout (uncomment to use)
# TOTAL_TIMEOUT=7200  # 2 hours total
# timeout $TOTAL_TIMEOUT bash -c '
# trap "echo \"Total timeout reached, killing all jobs\"; kill 0" EXIT

if [ $# -eq 0 ]; then
    echo "Usage: $0 /path/to/folder1/ [/path/to/folder2/] ..."
    echo "Example: $0 /media/disk_b/.../che_side0/"
    exit 1
fi

mkdir -p "$DESTINATION_DIR"
echo "Starting processing with up to $MAX_JOBS parallel jobs."
echo "Results will be copied to: $DESTINATION_DIR"
echo "========================================================"

all_subfolders=()
for FOLDER_PATH in "$@"; do
    FOLDER_PATH="${FOLDER_PATH%/}"
    if [ ! -d "$FOLDER_PATH" ]; then
        echo "Warning: Input folder '$FOLDER_PATH' does not exist. Skipping..."
        continue
    fi
    while IFS= read -r -d '' sub; do
        all_subfolders+=("$sub")
    done < <(find "$FOLDER_PATH" -mindepth 1 -maxdepth 2 -type d -print0)
done

echo "Found ${#all_subfolders[@]} subfolders to process."

job_count=0
for subfolder in "${all_subfolders[@]}"; do
    # IMPROVED job control
    while true; do
        # Clean up completed jobs first
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[i]}" 2>/dev/null; then
                unset "PIDS[i]"
            fi
        done
        # Reindex array to remove gaps
        PIDS=("${PIDS[@]}")
        
        if [ ${#PIDS[@]} -lt $MAX_JOBS ]; then
            break
        fi
        sleep 0.5
    done
    
    echo "[$(date '+%H:%M:%S')] Launching job $((job_count + 1))/${#all_subfolders[@]}: $(basename "$subfolder")"
    process_subfolder "$subfolder" "$DESTINATION_DIR" &
    PIDS+=($!)
    job_count=$((job_count + 1))
done

echo "========================================================"
echo "All $job_count jobs launched. Waiting for completion..."

# Better wait logic
while [ ${#PIDS[@]} -gt 0 ]; do
    for i in "${!PIDS[@]}"; do
        if ! kill -0 "${PIDS[i]}" 2>/dev/null; then
            wait "${PIDS[i]}" 2>/dev/null
            unset "PIDS[i]"
        fi
    done
    PIDS=("${PIDS[@]}")
    if [ ${#PIDS[@]} -gt 0 ]; then
        echo "[$(date '+%H:%M:%S')] Still waiting for ${#PIDS[@]} jobs to complete..."
        sleep 5
    fi
done

echo "All folders processed successfully!"
echo "Final results location: $DESTINATION_DIR"