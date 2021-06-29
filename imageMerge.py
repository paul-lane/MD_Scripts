from PIL import Image
# Based on code from:
# https://www.tutorialspoint.com/python_pillow/Python_pillow_merging_images.htm


#Read the images
image1 = Image.open('MDimages/C8mim0Side_crop.bmp')
image2 = Image.open('MDimages/C8mim10Side_crop.bmp')
image3 = Image.open('MDimages/C8mim25Side_crop.bmp')
image4 = Image.open('MDimages/C8mim50Side_crop.bmp')
image5 = Image.open('MDimages/C8mim75Side_crop.bmp')
image6 = Image.open('MDimages/C8mimFSide_crop.bmp')

# Assume all images are the same size as they were made with the image crop code
image1_size = image1.size

# Define new image to be 3* width of the original images and 2* the height
new_image = Image.new('RGB',(3*image1_size[0], 2*image1_size[1]), (250,250,250))

# Place image 1 in the top left corner
new_image.paste(image1,(0,0))
# Place image 2 next to image 1 (one image width away from top left)
new_image.paste(image2,(image1_size[0],0))
# Place image 3 next to image 2 (two image widths away from top left)
new_image.paste(image3,(2*image1_size[0],0))
# Place image 4 below image 1 (one image height away from top left)
new_image.paste(image4,(0,image1_size[1]))
# Place  image 5 below image 2
new_image.paste(image5,(image1_size[0],image1_size[1]))
# Place image 6 below image 3
new_image.paste(image6,(2*image1_size[0],image1_size[1]))

new_image.save("merged_image.bmp","BMP")
new_image.show()
