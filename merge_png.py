from PIL import Image
import os

def merge_images(image_dir, output_image_path, grid_size=(2, 5)):
    """
    Merges multiple images into a single image.
    
    Parameters:
    - image_dir: Directory containing the PNG images to be merged.
    - output_image_path: Path to save the output merged image.
    - grid_size: Tuple specifying the grid layout (rows, columns) for merging images.
    """
    # Get the list of all image file paths in the directory
    image_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.endswith('.png')]

    # Sort the image files for consistency (optional, based on your naming convention)
    image_files = sorted(image_files)

    # Select the first 10 images
    image_files = image_files[:10]

    # Open all images and find the maximum width and height
    images = [Image.open(img) for img in image_files]
    widths, heights = zip(*(img.size for img in images))

    max_width = max(widths)
    max_height = max(heights)

    # Define grid dimensions
    rows, cols = grid_size
    total_width = cols * max_width
    total_height = rows * max_height

    # Create a new blank image with the appropriate size
    merged_image = Image.new('RGB', (total_width, total_height))

    # Paste each image into the grid
    for index, img in enumerate(images):
        x_offset = (index % cols) * max_width
        y_offset = (index // cols) * max_height
        merged_image.paste(img, (x_offset, y_offset))

    # Save the merged image
    merged_image.save(output_image_path)
    print(f"Merged image saved at: {output_image_path}")

# Example usage:
image_directory = 'clade_plots'  # Directory where the 10 PNG images are saved
output_image_file = 'merged_clades_image.png'  # Output path for the merged image

# Merge images in a 2x5 grid (2 rows, 5 columns)
merge_images(image_directory, output_image_file, grid_size=(2, 5))
