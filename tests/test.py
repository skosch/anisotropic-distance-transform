import numpy as np
import matplotlib.pyplot as plt

import adt

# Create an image (here, just a square in the middle of zeros)
orig_img = np.zeros((200, 200)).astype(np.float32)
orig_img[90:110, 90:110] = 1.

# Define a maximum value (see the marching parabolas paper)
maxval = (np.sqrt(2)*200) ** 2 * 10 * 2

# Set everything to the maxval, except for the square (which becomes zero) ...
adt_img_1 = maxval - maxval * orig_img
# ... then apply the distance transform in place, like so
adt.adt(adt_img_1, 1, 1, 1, 1)

adt_img_2 = maxval - maxval * orig_img
# ... and, optionally, set factors for values for left, right, top, and bottom
adt.adt(adt_img_2, .5, 2., 1., 4.)

fig, ax = plt.subplots(1, 2)
ax[0].imshow(adt_img_1, cmap='gist_rainbow')
ax[0].contour(orig_img, cmap='gray')
ax[1].imshow(adt_img_2, cmap='gist_rainbow')
ax[1].contour(orig_img, cmap='gray')
plt.show()
