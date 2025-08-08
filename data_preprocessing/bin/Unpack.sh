#!/bin/bash

# Check if at least one folder path is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 /path/to/folder1/ [/path/to/folder2/] [/path/to/folder3/] [/path/to/folder4/]"
    echo "Example: $0 /media/disk_b/.../che_side0/ /media/disk_b/.../che_side1/"
    echo "You can provide 1-4 folders to process"
    exit 1
fi

# Function to process a single directory
process_directory() {
    local subfolder="$1"
    local subfolder_name=$(basename "$subfolder")
    
    # Check if calib.root file already exists - if so, skip this subfolder
    calib_exists=false
    for existing_calib in "${subfolder}"*calib.root; do
        if [ -e "$existing_calib" ]; then
            calib_exists=true
            break
        fi
    done
    
    if [ "$calib_exists" = true ]; then
        echo "⏭ Skipping subfolder: $subfolder_name (calib.root already exists)"
        return
    fi
    
    echo "Processing subfolder: $subfolder_name"
    
    # Count .daq files first
    daq_count=0
    for f in "${subfolder}"*.daq; do
        if [ -e "$f" ]; then
            ((daq_count++))
        fi
    done
    
    # First, run TDMunpack on all .daq files in the subfolder (only if 4 or fewer)
    if [ $daq_count -gt 4 ]; then
        echo "  ⚠ Found $daq_count .daq files (>4) - Skipping TDMunpack step"
    else
        echo "  Running TDMunpack on .daq files in $subfolder_name ($daq_count files)"
        daq_files_found=false
        for f in "${subfolder}"*.daq; do
            # Check if the glob matched any files
            if [ -e "$f" ]; then
                daq_files_found=true
                echo "    Processing: $(basename "$f")"
                if ./TDMunpack -f "$f"; then
                    echo "    ✓ TDMunpack successful for $(basename "$f")"
                else
                    echo "    ✗ TDMunpack failed for $(basename "$f") (exit code: $?)"
                fi
            fi
        done
        
        if [ "$daq_files_found" = false ]; then
            echo "  → No .daq files found in $subfolder_name"
        fi
    fi
    
    # Run the ls command to create febs_files_list.list
    echo "  Running: ls ${subfolder}*Slot* > febs_files_list.list"
    if ls "${subfolder}"*Slot* > febs_files_list.list 2>/dev/null; then
        echo "  ✓ File list created successfully"
        
        # Run the unpack command
        echo "  Running: ./unpack"
        if ./unpack; then
            echo "  ✓ Unpack completed successfully for $subfolder_name"
        else
            echo "  ✗ Unpack failed for $subfolder_name (exit code: $?)"
            echo "  → Continuing to next subfolder..."
        fi
    else
        echo "  ✗ No *Slot* files found in $subfolder_name or ls command failed"
        echo "  → Skipping ./unpack for this subfolder"
    fi
    
    raw_files_found=false
    for f in "${subfolder}"*raw.root; do
        # Check if the glob matched any files
        if [ -e "$f" ]; then
            raw_files_found=true
            echo "    Processing: $(basename "$f")"
            if ./Calibration "$f" "55" "50"; then
                echo "    ✓ Calib successful for $(basename "$f")"
            else
                echo "    ✗ Calib failed for $(basename "$f") (exit code: $?)"
            fi
        fi
    done
    
    if [ "$raw_files_found" = false ]; then
        echo "  → No raw.root files found in $subfolder_name"
    fi

    calib_files_found=false
    for f in "${subfolder}"*calib.root; do
        # Check if the glob matched any files
        if [ -e "$f" ]; then
            calib_files_found=true
            echo "    Processing: $(basename "$f")"
            if ./EventStructure_12free "$f"; then
                echo "    ✓ EventStructure successful for $(basename "$f")"
            else
                echo "    ✗ EventStructure failed for $(basename "$f") (exit code: $?)"
            fi
        fi
    done
    
    if [ "$calib_files_found" = false ]; then
        echo "  → No calib.root files found in $subfolder_name"
    fi

    echo "  ────────────────────────────────"
}

# Process each provided folder
for FOLDER_PATH in "$@"; do
    # Remove trailing slash if present
    FOLDER_PATH="${FOLDER_PATH%/}"
    
    # Check if the folder exists
    if [ ! -d "$FOLDER_PATH" ]; then
        echo "Error: Folder '$FOLDER_PATH' does not exist. Skipping..."
        continue
    fi
    
    echo ""
    echo "========================================"
    echo "Processing subfolders in: $FOLDER_PATH"
    echo "========================================"
    
    # First, try level 1 subdirectories and check if any have processable files
    level1_dirs_with_files=0
    level1_dirs=()
    
    # Collect level 1 directories and check for files
    for subfolder in "$FOLDER_PATH"/*/; do
        if [ -d "$subfolder" ]; then
            level1_dirs+=("$subfolder")
            
            # Check if this directory has any relevant files (.daq, *Slot*, *raw.root, *calib.root)
            has_files=false
            for pattern in "*.daq" "*Slot*" "*raw.root" "*calib.root"; do
                for file in "${subfolder}"$pattern; do
                    if [ -e "$file" ]; then
                        has_files=true
                        break 2
                    fi
                done
            done
            
            if [ "$has_files" = true ]; then
                ((level1_dirs_with_files++))
            fi
        fi
    done
    
    # If no level 1 directories have processable files, go to level 2
    if [ $level1_dirs_with_files -eq 0 ] && [ ${#level1_dirs[@]} -gt 0 ]; then
        echo "No processable files found in level 1 directories. Searching level 2..."
        
        # Process level 2 subdirectories
        for level1_dir in "${level1_dirs[@]}"; do
            for subfolder in "$level1_dir"/*/; do
                if [ -d "$subfolder" ]; then
                    process_directory "$subfolder"
                fi
            done
        done
    else
        # Process level 1 directories
        if [ ${#level1_dirs[@]} -eq 0 ]; then
            echo "No subdirectories found in $FOLDER_PATH"
        else
            echo "Found $level1_dirs_with_files level 1 directories with processable files. Processing level 1..."
            for subfolder in "${level1_dirs[@]}"; do
                process_directory "$subfolder"
            done
        fi
    fi
    
    echo "Completed processing: $FOLDER_PATH"
done

echo ""
echo "All folders processed."