from pathlib import Path
import argparse
from PIL import Image


IMAGE_EXTENSIONS = {".jpg", ".jpeg", ".png", ".bmp", ".gif", ".tiff", ".webp"}


def combine_images(first_path, second_path, output_path):
    first = Image.open(first_path).convert("RGB")
    second = Image.open(second_path).convert("RGB")

    # Make both images the same width by resizing the second image.
    if first.width != second.width:
        new_height = round(second.height * first.width / second.width)
        second = second.resize((first.width, new_height), Image.Resampling.LANCZOS)

    # Create a canvas tall enough for both images.
    combined = Image.new(
        "RGB",
        (first.width, first.height + second.height),
        "white",
    )

    combined.paste(first, (0, 0))
    combined.paste(second, (0, first.height))

    combined.save(output_path)

def create_combined_images(first_folder, second_folder, output_folder):

    first_folder = Path(first_folder)
    second_folder = Path(second_folder)
    output_folder = Path(output_folder)


    output_folder.mkdir(parents=True, exist_ok=True)

    first_images = {
        file.name: file
        for file in first_folder.iterdir()
        if file.is_file() and file.suffix.lower() in IMAGE_EXTENSIONS
    }

    second_images = {
        file.name: file
        for file in second_folder.iterdir()
        if file.is_file() and file.suffix.lower() in IMAGE_EXTENSIONS
    }

    matching_names = sorted(first_images.keys() & second_images.keys())

    print(f"Found {len(matching_names)} matching image(s).")

    for name in matching_names:
        output_file = output_folder / name

        print(f"Combining: {name}")

        try:
            combine_images(
                first_images[name],
                second_images[name],
                output_file,
            )
        except Exception as e:
            print(f"  Error processing {name}: {e}")

    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description="Combine matching images vertically from two folders."
    )
    parser.add_argument("first_folder", type=Path)
    parser.add_argument("second_folder", type=Path)
    parser.add_argument("output_folder", type=Path)

    args = parser.parse_args()

    create_combined_images(args.first_folder, args.second_folder, args.output_folder)

if __name__ == "__main__":
    main()