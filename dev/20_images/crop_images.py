from pathlib import Path
import cv2
import numpy as np


def find_non_whitespace_pixels(image):
    # Convert the image to grayscale to check for white pixels
    grayscale_image = np.sum(np.abs(image - [255, 255, 255]), axis=-1)

    # Find the first non-white pixel from the left
    left = np.argmax(np.sum(grayscale_image > 0, axis=1) > 0)

    # Find the first non-white pixel from the right
    right = grayscale_image.shape[0] - np.argmax(np.sum(grayscale_image[::-1] > 0, axis=1) > 0)

    # Find the first non-white pixel from the top
    top = np.argmax(np.sum(grayscale_image > 0, axis=0) > 0)

    # Find the first non-white pixel from the bottom
    bottom = grayscale_image.shape[1] - np.argmax(np.sum(grayscale_image[:, ::-1] > 0, axis=0) > 0)

    return left, right, top, bottom


if __name__ == "__main__":
    image_dir = Path(Path.home(), "Documents/xray_ms/figure_3/7mhf_30")
    crop_image_dir = Path(Path.home(), "Documents/xray_ms/figure_3/{}_0".format(image_dir.name))
    crop_image_dir.mkdir(exist_ok=True)

    for image_file in image_dir.glob("*.png"):
        print(image_file)
        img = cv2.imread(str(image_file))

        left, right, top, bottom = find_non_whitespace_pixels(img)

        print(left, right, top, bottom)

        crop_img = img[left:right+1,top:bottom+1,]
        cv2.imwrite(str(Path(crop_image_dir, image_file.name)), crop_img)
