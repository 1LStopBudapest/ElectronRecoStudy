import os
from PIL import Image, ImageChops
import re

def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('f1', type=str)
    argParser.add_argument('f2', type=str)
    argParser.add_argument('f3', type=str)
    return argParser

def crop_white_borders(image, bg_color=(255, 255, 255), border=5):
    """Crop white borders from the image but leave a small white border (default: 5 pixels)."""
    bg = Image.new(image.mode, image.size, bg_color)
    diff = ImageChops.difference(image, bg)
    bbox = diff.getbbox()
    
    if bbox:
        left = max(bbox[0] - border, 0)
        upper = max(bbox[1] - border, 0)
        right = min(bbox[2] + border, image.width)
        lower = min(bbox[3] + border, image.height)
        return image.crop((left, upper, right, lower))
    
    return image  # if no border found

def load_and_crop_images(folder_path, image_names):
    """Load and crop the images from a folder."""
    images = []
    for name in image_names:
        path = os.path.join(folder_path, name)
        img = Image.open(path).convert("RGB")
        img = crop_white_borders(img)
        images.append(img)
    return images

def create_2x2_grid(images):
    """Arrange 4 images in a 2x2 grid."""
    assert len(images) == 4, "Need exactly 4 images"
    
    # Resize all to same size (e.g., min width/height among them)
    min_width = min(img.width for img in images)
    min_height = min(img.height for img in images)
    resized = [img.resize((min_width, min_height)) for img in images]

    grid_width = min_width * 2
    grid_height = min_height * 2
    grid = Image.new('RGB', (grid_width, grid_height), (255, 255, 255))

    # Paste images
    grid.paste(resized[0], (0, 0))
    grid.paste(resized[1], (min_width, 0))
    grid.paste(resized[2], (0, min_height))
    grid.paste(resized[3], (min_width, min_height))

    return grid

# === Main setup ===


options = get_parser().parse_args()

folder1 = options.f1
folder2 = options.f2


def natural_sort_key(s):
    """Helper function to sort strings with numbers naturally."""
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]  
pattern = re.compile(r'^figure(\d+)\.png$', re.IGNORECASE)

png_files = sorted(
    [f for f in os.listdir(options.f1) if pattern.match(f)],
    key=natural_sort_key
)


for i in range(0,len(png_files),2):

    # Image names (must exist in both folders)
    image_names = [png_files[i], png_files[i+1]]

    # Load and crop
    images_folder1 = load_and_crop_images(folder1, image_names)
    images_folder2 = load_and_crop_images(folder2, image_names)

    # Combine all images (2 from each folder)
    all_images = images_folder1 + images_folder2

    # Create and save the grid
    grid_image = create_2x2_grid(all_images)
    grid_image.save(options.f3+"/figure"+str(i)+".png")
    #grid_image.show()

