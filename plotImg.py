import cv2
import matplotlib.pyplot as plt

# Function to display the image and capture the clicked pixel value
def on_click(event, x, y, flags, param):
    if event == cv2.EVENT_LBUTTONDOWN:  # Left mouse button clicked
        # Print the clicked pixel coordinates and its value
        pixel_value = img[y, x]
        print(f"Clicked at ({x}, {y}) - Pixel Value: {pixel_value}")

# Load the image
img = cv2.imread('scenes/cbox/cbox-marble.png')

# Create a window and set the mouse callback function
cv2.namedWindow('Image')
cv2.setMouseCallback('Image', on_click)

# Display the image
while True:
    cv2.imshow('Image', img)
    if cv2.waitKey(20) & 0xFF == 27:  # Press 'Esc' to exit
        break

# Destroy all windows when done
cv2.destroyAllWindows()
