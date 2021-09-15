# based on: 
# https://www.geeksforgeeks.org/python-pil-image-crop-method/

# Importing Image class from PIL module
from PIL import Image

# Opens a image in RGB mode
im = Image.open(r"MDimages/C8mimFSide.bmp")
 
# Size of the image in pixels (size of original image)
width, height = im.size
 
# Setting the points for cropped image
left = 5 * width / 16
top = 3 * height / 16
right = 11 * width /16
bottom = 13 * height / 16
 
# Cropped image of above dimension
# (It will not change original image)
im1 = im.crop((left, top, right, bottom))
 
# Shows the image in image viewer
im1.show()

# Save image
im1 = im1.save("MDimages/C8mimFSide_crop.bmp") 
