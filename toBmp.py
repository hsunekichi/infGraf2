import os
from PIL import Image

def convert_to_png(input_dir, output_dir):
    """
    Converts all images in the input directory to PNG format and saves them in the output directory.

    Args:
        input_dir (str): Path to the input directory containing images.
        output_dir (str): Path to the output directory to save converted PNG images.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        input_path = os.path.join(input_dir, filename)
        
        # Skip directories
        if os.path.isdir(input_path):
            continue

        # Open and convert the file if it's an image
        try:
            with Image.open(input_path) as img:
                # Convert the image filename to .bmp
                base_name = os.path.splitext(filename)[0]
                output_path = os.path.join(output_dir, f"{base_name}.bmp")
                
                # Save as PNG
                img.save(output_path, format="BMP")
                print(f"Converted: {input_path} -> {output_path}")
        except Exception as e:
            print(f"Skipping {input_path}: {e}")

if __name__ == "__main__":
    input_dir = input("Enter the path to the input directory: ")
    output_dir = input("Enter the path to the output directory: ")
    convert_to_png(input_dir, output_dir)
