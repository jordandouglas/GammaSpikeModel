#!/bin/bash

# --- R Script Executor (Centralized) ---
# This script runs a central R script across multiple analysis folders.

# Exit immediately if a command exits with a non-zero status.
set -e

# --- 1. DEFINE CENTRAL SCRIPT PATH ---
# The location of our single, authoritative R script.
# This path is relative to where this shell script is run.
SCRIPT_DIR=$(dirname "$0")
R_SCRIPT_PATH="$SCRIPT_DIR/summarize_posteriors.R"

# --- 2. PARSE & VALIDATE ARGUMENTS ---
# (This entire section remains unchanged from your original script)
BURNIN=""
NAME_TAG=""
FOLDERS=()
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --burnin)
            if [ -n "$2" ]; then BURNIN="$2"; shift 2; else echo "❌ Error: --burnin requires a value." >&2; exit 1; fi ;;
        --name)
            if [ -n "$2" ]; then NAME_TAG="$2"; shift 2; else echo "❌ Error: --name requires a value." >&2; exit 1; fi ;;
        *)
            FOLDERS+=("$1"); shift ;;
    esac
done

USAGE_MSG="Usage: ./run_posterior_summaries.sh --burnin <value> [--name <tag>] <folder1> [folder2] ..."
if [ -z "$BURNIN" ] || [ ${#FOLDERS[@]} -eq 0 ]; then
    echo "❌ Error: Missing required arguments."
    echo "$USAGE_MSG"
    exit 1
fi

# --- 3. PROCESS EACH FOLDER ---
echo "🚀 Starting batch processing with burn-in = $BURNIN..."

for FOLDER in "${FOLDERS[@]}"; do
    echo "----------------------------------------"
    echo "➡️ Processing folder: '$FOLDER'"

    # --- KEY CHANGE: Use the FOLDER as the working directory ---
    if [ ! -d "$FOLDER" ]; then
        echo "⚠️  Warning: Folder '$FOLDER' not found. Skipping."
        continue
    fi

    TEMPLATES_DIR="$FOLDER/templates"
    if [ ! -d "$TEMPLATES_DIR" ]; then
        echo "⚠️  Warning: Templates directory not found in '$FOLDER'. Skipping."
        continue
    fi
    
    # Calculate Parameters
    NSIMS=$(find "$TEMPLATES_DIR" -mindepth 1 -maxdepth 1 -type d -printf '.' | wc -c)
    echo "   Found $NSIMS simulation subfolders."

    # Construct the output file name
    FOLDER_BASENAME=$(basename "$FOLDER")
    if [ -n "$NAME_TAG" ]; then
        OUTPUT_FILE="posterior_summary_${NAME_TAG}_${FOLDER_BASENAME}.rds"
    else
        OUTPUT_FILE="posterior_summary_${FOLDER_BASENAME}.rds"
    fi
    echo "   Output file will be '$OUTPUT_FILE'."

    # --- KEY CHANGE: Execute from within the folder ---
    # We use a subshell `()` so the 'cd' command only affects this block.
    # The R script is called using its central path.
    # The output path no longer needs '../' because we are in the right place.
    (
        echo "   Changing working directory to '$FOLDER'..."
        cd "$FOLDER"
        
        Rscript "../$R_SCRIPT_PATH" \
            --nsims=$NSIMS \
            --burnin=$BURNIN \
            --out="$OUTPUT_FILE"
    )

    echo "✅ Finished processing '$FOLDER'."
done

echo "----------------------------------------"
echo "🎉 All tasks completed!"