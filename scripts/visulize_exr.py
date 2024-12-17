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


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python visualize_exr.py <path_to_exr_file>")
        sys.exit(1)

    image_path = sys.argv[1]
    if not os.path.exists(image_path):
        print(f"File not found: {image_path}")
        sys.exit(1)

    display_exr_with_reinhard(image_path)