
import cv2
import numpy as np

# Load the image
image_path = '/Users/zongmin/Downloads/mask5.jpg'
image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

# Binarize the image (assume white is the structure and black is the background)
_, binary_image = cv2.threshold(image, 127, 255, cv2.THRESH_BINARY)

# Get the total width and height of the image
total_height = binary_image.shape[0]
total_width = binary_image.shape[1]

# Number of columns to divide the image into
num_columns = 20
column_width = total_width // num_columns

# Initialize the SWC content
swc_content = []
node_id = 1
parent_id = -1
y_center = total_height // 2  # Fixed y position

# Scan each of the 20 columns
for i in range(num_columns):
    x_start = i * column_width
    x_end = (i + 1) * column_width if i < num_columns - 1 else total_width

    # Find the white regions in this column segment (all white pixels)
    y_indices = np.where(binary_image[:, x_start:x_end].any(axis=1))[0]
    if len(y_indices) > 0:
        # Calculate the maximum connected white region's length
        y_min = np.min(y_indices)
        y_max = np.max(y_indices)
        z_center = (y_min + y_max) // 2  # Z coordinate is the center of the white region
        diameter = y_max - y_min  # The length of the white region (to be used as diameter)

        # Set the x coordinate to the middle of the current column
        x_center = (x_start + x_end) // 2

        # Add to SWC content: (node_id, type, x_center, y_center, z_center, diameter/2, parent_id)
        swc_content.append([node_id, 2, x_center, y_center, z_center, diameter / 2.0, parent_id])
        parent_id = node_id
        node_id += 1

# Save the SWC content to a file
swc_output_path = '/Users/zongmin/Downloads/1.swc'
with open(swc_output_path, 'w') as f:
    f.write('# SWC file generated from image\n')
    f.write('# inode R X Y Z D/2 idpar\n')
    for line in swc_content:
        f.write('\t'.join(map(str, line)) + '\n')

print(f"SWC file saved to {swc_output_path}")
