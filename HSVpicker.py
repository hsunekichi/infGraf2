import cv2
import numpy as np

# Function to be called when the mouse is clicked
def click_event(event, x, y, flags, param):
    if event == cv2.EVENT_LBUTTONDOWN:  # On left mouse button click
        # Get the BGR value of the pixel
        bgr_pixel = image[y, x]
        
        # Convert the BGR pixel to HSV
        hsv_pixel = cv2.cvtColor(np.uint8([[bgr_pixel]]), cv2.COLOR_BGR2HSV)[0][0]
        
        # Print the HSV values
        print(f'HSV value at ({x}, {y}): {hsv_pixel}')

# Load the image
image_path = 'LedaOriginal.jpg'
image = cv2.imread(image_path)

if image is None:
    print("Error: Could not open or find the image.")
else:
    # Create a window and set the mouse callback function
    cv2.imshow('Image', image)
    cv2.setMouseCallback('Image', click_event)
    
    # Wait until a key is pressed and destroy all windows
    cv2.waitKey(0)
    cv2.destroyAllWindows()
