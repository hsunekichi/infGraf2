import OpenEXR
import Imath
import sys
import numpy as np

# Function to read an EXR file and return the pixel data as a NumPy array
def read_exr(file_path):
    exr_file = OpenEXR.InputFile(file_path)
    
    # Get header information
    header = exr_file.header()
    
    # Get the data window
    dw = header['dataWindow']
    width = dw.max.x - dw.min.x + 1
    height = dw.max.y - dw.min.y + 1
    
    # Define the channel format (assuming RGB, and FLOAT channels)
    FLOAT = Imath.PixelType(Imath.PixelType.FLOAT)
    
    # Read each channel (R, G, B) into separate arrays
    red_str = exr_file.channel('R', FLOAT)
    green_str = exr_file.channel('G', FLOAT)
    blue_str = exr_file.channel('B', FLOAT)
    
    # Convert strings to numpy arrays
    red = np.frombuffer(red_str, dtype=np.float32).reshape(height, width)
    green = np.frombuffer(green_str, dtype=np.float32).reshape(height, width)
    blue = np.frombuffer(blue_str, dtype=np.float32).reshape(height, width)
    
    # Combine the channels into a single array
    img = np.stack([red, green, blue], axis=-1)
    
    return img

# Function to compute RMSE between two images
def compute_rmse(img_ref, img_pred):
    # Ensure both images have the same shape
    if img_ref.shape != img_pred.shape:
        raise ValueError("Images must have the same dimensions")
    
    # Compute the RMSE
    mse = np.mean((img_ref - img_pred) ** 2)
    rmse = np.sqrt(mse)
    
    return rmse


if len(sys.argv) < 2:
    print("Usage: python exrMetrics.py <exr_file_predict>")
    sys.exit(1)

# Example usage
exr_file_ref = "scenes/pa4/cbox/cbox-whittedRef.exr"
exr_file2 = sys.argv[1]

# Load EXR images
img_ref = read_exr(exr_file_ref)
img2 = read_exr(exr_file2)

# Compute RMSE
rmse_value = compute_rmse(img_ref, img2)
print(f"RMSE between the images: {rmse_value}")
