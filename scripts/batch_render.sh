#!/bin/bash


if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <scene_path> <output_dir> [num_iterations]"
    echo "  scene_path: Path to the input scene XML"
    echo "  output_dir: Directory to save output EXR files"
    echo "  num_iterations: Number of iterations"
    exit 1
fi

# Configuration
SCENE_PATH="$1"      # Path to the input scene XML
OUT_DIR="$2"         # Directory to save output EXR files
NUM_ITERATIONS=${3:-1}  # Number of iterations

# Ensure OUT_DIR exists
mkdir -p "$OUT_DIR"

# Get highest output iteration number
if [ -z "$(ls -A $OUT_DIR)" ]; then
    i=1
else
    i=$(ls -1 $OUT_DIR | grep -o '[0-9]*' | sort -n | tail -1)
    i=$((i+1))
fi

# Run the render engine multiple times
# While i == -1 or i <= NUM_ITERATIONS
while [[  "$#" -eq 2  || $i -le $NUM_ITERATIONS ]]; do
    echo "Starting iteration $i..."
    
    # Run the render engine
    echo "Rendering scene: $SCENE_PATH"
    ./build/nori "$SCENE_PATH" --nogui
    
    # Construct input and output file paths
    EXR_FILE="${SCENE_PATH%.xml}.exr"             # Output EXR path
    EXR_OUTPUT="${OUT_DIR}/output_${i}.exr"       # Output file with iteration ID
    
    # Verify the render engine created the EXR file
    if [ -f "$EXR_FILE" ]; then
        # Copy the EXR file to the output directory with the iteration ID
        cp "$EXR_FILE" "$EXR_OUTPUT"
        echo "Saved iteration $i to: $EXR_OUTPUT"
        
        # Optional: Remove the original EXR file
        rm "$EXR_FILE"
    else
        echo "Error: Expected EXR file not found after iteration $i."
        exit 1
    fi

    i=$((i+1))
done

echo "Batch render complete. Files saved to: $OUT_DIR"
