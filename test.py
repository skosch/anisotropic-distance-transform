
import numpy as np

import C_gdt

import matplotlib.pyplot as plt

# Create a random image

orig_img = np.zeros((100, 200)).astype(np.float32)
maxval = (np.sqrt(2)*200) ** 2 * 3 * 2
orig_img[40:60, 40:60] = 1.

gdt_img = maxval - maxval * orig_img
C_gdt.gdt(gdt_img, 2., 3., 4, 5.)
#C_gdt.gdt(gdt_img, 1, 1, 1, 1)
plt.imshow(orig_img, cmap='gray')
plt.imshow(gdt_img, alpha=0.2, cmap='gray')
plt.imshow(gdt_img, alpha=0.2, cmap='inferno')
plt.contour(gdt_img, 20)
plt.show()
