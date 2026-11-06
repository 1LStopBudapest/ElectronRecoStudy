import os
import argparse
from PIL import Image

def create_combined_images(f1, f2, f3):

    os.makedirs(f3, exist_ok=True)

    # Get and sort f2 images starting with binX
    f2_images = [f for f in os.listdir(f2) if f.startswith("bin") and f.endswith(".png")]
    f2_images.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))

    # Get and sort f1 images starting with binX
    f1_images = [f for f in os.listdir(f1) if f.startswith("bin") and f.endswith(".png")]
    f1_images.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))

    # Pair images
    for idx, f2_img_name in enumerate(f2_images):
        if idx >= len(f1_images):
            print("Warning: Not enough f1 images for "+str(f2_img_name))
            break
        
        f2_path = os.path.join(f2, f2_img_name)
        f1_path = os.path.join(f1, f1_images[idx])

        img1 = Image.open(f2_path)
        img2 = Image.open(f1_path)

        # Determine new image size (width=max width, height=sum of heights)
        new_width = max(img1.width, img2.width)
        new_height = img1.height + img2.height

        # Create a white background image
        new_img = Image.new("RGB", (new_width, new_height), (255, 255, 255))

        # Paste img2 centered on top
        top_x = (new_width - img2.width) // 2
        new_img.paste(img2, (top_x, 0))

        # Paste img1 centered on bottom
        bottom_x = (new_width - img1.width) // 2
        new_img.paste(img1, (bottom_x, img2.height))

        # Save to f3 with the f2 filename
        save_path = os.path.join(f3, f2_img_name)
        new_img.save(save_path)

    print("All images processed! Output: "+str(f3))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine images from two folders into a third folder.")
    parser.add_argument("f1", help="Path to first folder (f1) containing template DATA")
    parser.add_argument("f2", help="Path to second folder (f2) containing MC")
    parser.add_argument("f3", help="Path to output folder (f3)")

    args = parser.parse_args()
    create_combined_images(args.f1, args.f2, args.f3)