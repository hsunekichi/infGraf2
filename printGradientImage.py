import OpenEXR
import Imath
import numpy as np
import cv2

# Function to read the EXR image and return the RGB channels
def load_exr_image(file_path):
    exr_file = OpenEXR.InputFile(file_path)
    header = exr_file.header()
    
    dw = header['dataWindow']
    width = dw.max.x - dw.min.x + 1
    height = dw.max.y - dw.min.y + 1
    
    # Read in the RGB channels
    FLOAT = Imath.PixelType(Imath.PixelType.FLOAT)
    r_str = exr_file.channel('R', FLOAT)
    g_str = exr_file.channel('G', FLOAT)
    b_str = exr_file.channel('B', FLOAT)
    
    # Convert the strings to numpy arrays
    R = np.frombuffer(r_str, dtype=np.float32).reshape((height, width))
    G = np.frombuffer(g_str, dtype=np.float32).reshape((height, width))
    B = np.frombuffer(b_str, dtype=np.float32).reshape((height, width))
    
    return R, G, B

# Function to compute luminance from RGB
def compute_luminance(R, G, B):
    # Luminance formula for sRGB space
    luminance = 0.2126 * R + 0.7152 * G + 0.0722 * B
    return luminance

# Function to compute gradients using Sobel operator
def compute_gradients(luminance):
    # Sobel operator to compute x and y gradients
    grad_x = cv2.Sobel(luminance, cv2.CV_64F, 1, 0, ksize=3)  # Gradient in x
    grad_y = cv2.Sobel(luminance, cv2.CV_64F, 0, 1, ksize=3)  # Gradient in y

    # Compute gradient magnitude
    gradient_magnitude = np.sqrt(grad_x**2 + grad_y**2)
    
    return grad_x, grad_y, gradient_magnitude

def compute_gradients_finite_difference(luminance):
    # Compute gradients using finite differences
    grad_x = np.gradient(luminance, axis=1)
    grad_y = np.gradient(luminance, axis=0)

    # Compute gradient magnitude
    gradient_magnitude = np.sqrt(grad_x**2 + grad_y**2)
    
    return grad_x, grad_y, gradient_magnitude

# Main function to load EXR, compute luminance, and compute gradients
def process_exr(file_path):
    # Step 1: Load the EXR image
    R, G, B = load_exr_image(file_path)
    
    # Step 2: Compute the luminance matrix
    luminance = compute_luminance(R, G, B)
    
    # Step 3: Compute the gradients
    grad_x, grad_y, gradient_magnitude = compute_gradients(luminance)
    
    return luminance, grad_x, grad_y, gradient_magnitude

if __name__ == "__main__":
    # Path to your EXR image
    exr_file_path = 'scenes/pa4/cbox/cbox-whitted.exr'

    # Process the EXR image
    luminance, grad_x, grad_y, gradient_magnitude = process_exr(exr_file_path)

    # Print or visualize results
    # Plot image in greyscale
    cv2.imshow('Sobel gradient magnitude', gradient_magnitude)
    
    cv2.waitKey(0)