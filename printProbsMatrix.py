import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Usage: python printProbsMatrix.py <numberSPP>")
    sys.exit(1)

numberSPP = sys.argv[1]

# Load the matrix from the text file
matrix = np.loadtxt(f'spp{numberSPP}.txt').T

# Print the loaded matrix
print("Matrix heigth: ", len(matrix), "Matrix width: ", len(matrix[0]))

# Plotting the matrix as a black-and-white (grayscale) image
plt.imshow(matrix, cmap='gray', interpolation='nearest')

# Optional: Remove axis ticks for cleaner visualization
plt.xticks([])
plt.yticks([])

# Show the image
plt.show()