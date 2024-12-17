import os
import numpy as np
import OpenEXR
import Imath
import imageio.v2 as imageio
from PIL import Image
import cv2
import sys

def load_exr(filepath):
    """Load an EXR file and return its pixel data as a numpy array."""
    exr_file = OpenEXR.InputFile(filepath)
    header = exr_file.header()
    
    # Get image size
    dw = header['dataWindow']
    width = dw.max.x - dw.min.x + 1
    height = dw.max.y - dw.min.y + 1

    # Extract channel names
    channel_names = header['channels'].keys()  # Get channel names from the header
    
    # Select only the first three channels (RGB) if available
    rgb_channels = ['R', 'G', 'B']
    available_channels = [c for c in rgb_channels if c in channel_names]
    if len(available_channels) < 3:
        raise ValueError(f"Image at {filepath} does not contain all RGB channels.")

    # Read all RGB channels and convert to numpy array
    channel_data = [np.frombuffer(exr_file.channel(c, Imath.PixelType(Imath.PixelType.FLOAT)), 
                                  dtype=np.float32).reshape((height, width)) for c in available_channels]
    exr_file.close()
    
    # Stack RGB channels
    return np.stack(channel_data, axis=-1)  # Shape: (H, W, 3)


def save_exr(filepath, data):
    """Save pixel data to an EXR file."""
    height, width, _ = data.shape
    header = OpenEXR.Header(width, height)
    channels = {'R': data[:, :, 0].astype(np.float32).tobytes(),
                'G': data[:, :, 1].astype(np.float32).tobytes(),
                'B': data[:, :, 2].astype(np.float32).tobytes()}
    
    exr_file = OpenEXR.OutputFile(filepath, header)
    exr_file.writePixels(channels)
    exr_file.close()


def display_exr_with_reinhard(image_path):
    """Display an EXR image with interactive Reinhard tone mapping sliders."""
    # Load the HDR image
    hdr_image = load_exr(image_path)

    # Window for display
    cv2.namedWindow("Reinhard Tone Mapping", cv2.WINDOW_NORMAL)

    # Reinhard parameters (initial values)
    params = {"gamma": 24, "intensity": 0, "light_adapt": 70, "color_adapt": 70}

    # Callback function for trackbar (does nothing, needed by OpenCV)
    def on_trackbar(val):
        pass

    # Create trackbars for Reinhard tone mapping parameters
    cv2.createTrackbar("Gamma x10", "Reinhard Tone Mapping", params["gamma"], 50, on_trackbar)
    cv2.createTrackbar("Intensity", "Reinhard Tone Mapping", params["intensity"], 100, on_trackbar)
    cv2.createTrackbar("Light Adapt", "Reinhard Tone Mapping", params["light_adapt"], 100, on_trackbar)
    cv2.createTrackbar("Color Adapt", "Reinhard Tone Mapping", params["color_adapt"], 100, on_trackbar)

    while True:
        # Read trackbar positions
        gamma = cv2.getTrackbarPos("Gamma x10", "Reinhard Tone Mapping") / 10.0
        intensity = (cv2.getTrackbarPos("Intensity", "Reinhard Tone Mapping") / 100.0) * 16 - 8 # Range -8, 8
        light_adapt = cv2.getTrackbarPos("Light Adapt", "Reinhard Tone Mapping") / 100.0
        color_adapt = cv2.getTrackbarPos("Color Adapt", "Reinhard Tone Mapping") / 100.0

        # Create Reinhard tone mapper with current parameters
        tonemap = cv2.createTonemapReinhard(gamma=gamma, intensity=intensity,
                                            light_adapt=light_adapt, color_adapt=color_adapt)
        ldr_image = tonemap.process(hdr_image)

        # Clip and convert the result to 8-bit for display
        ldr_image = np.clip(ldr_image * 255, 0, 255).astype(np.uint8)

        # Invert channels
        ldr_image = cv2.cvtColor(ldr_image, cv2.COLOR_BGR2RGB)

        # Show the image
        cv2.imshow("Reinhard Tone Mapping", ldr_image)

        # Break on ESC key
        if cv2.waitKey(1) & 0xFF == 27:
            break

    cv2.destroyAllWindows()

def gamma_tone_mapping(data, gamma=2.4):
    """Apply gamma correction to the pixel data."""
    return np.clip(data ** (1.0 / gamma), 0, 1)



def reinhard_tone_mapping(hdr_image):
    """Apply Reinhard global tone mapping to HDR image."""
    tonemap = cv2.createTonemapReinhard(gamma=2.2, intensity=0, light_adapt=0.7, color_adapt=0.7)
    ldr = tonemap.process(hdr_image)  # Perform tone mapping
    return np.clip(ldr, 0, 1)  # Ensure values are in [0, 1]




def main(input_dir, output_exr, output_png):
    """Main function to combine EXR images and export them."""
    # Get list of EXR files
    exr_files = [f for f in os.listdir(input_dir) if f.endswith('.exr')]
    if not exr_files:
        print("No EXR files found in the directory.")
        return

    print(f"Found {len(exr_files)} EXR files. Loading...")
    
    # Load and sum all EXR images
    accumulated = None
    for idx, file in enumerate(exr_files):
        filepath = os.path.join(input_dir, file)
        data = load_exr(filepath)
        print(f"Loaded: {file}")
        if accumulated is None:
            accumulated = data
        else:
            accumulated += data

    # Compute the mean
    mean_image = accumulated / len(exr_files)
    print("Computed the mean of all images.")

    # Save the EXR file
    save_exr(output_exr, mean_image)
    print(f"Saved mean image to: {output_exr}")

    display_exr_with_reinhard(output_exr)

    # Apply gamma correction and save as PNG
    gamma_corrected = reinhard_tone_mapping(mean_image)
    gamma_corrected_8bit = (gamma_corrected * 255).astype(np.uint8)
    Image.fromarray(gamma_corrected_8bit).save(output_png)
    print(f"Saved gamma-corrected image to: {output_png}")

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python combineEXR.py <input_directory> <output_file_prefix>")
        sys.exit(1)

    # Edit these paths as needed
    input_directory = sys.argv[1]  # Directory containing EXR files
    output_exr_path = sys.argv[2]+".exr"
    output_png_path = sys.argv[2]+".png"

    main(input_directory, output_exr_path, output_png_path)
