# Image Compressor
# Adapted from Brunton, Kutz - databookV2

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
import matplotlib.cm as cm

A = imread('ahmed.jpg')
X = np.mean(A, -1)

plot = plt.imshow(X, cmap = cm.gray)

U, S, VT = np.linalg.svd(X, full_matrices = False)
S = np.diag(S)

for r in (2000, 500, 150):
    alpha = int(np.floor(np.min(np.shape(X)) / r))
    Xapprox = U[:,:alpha] @ S[0:alpha,:alpha] @ VT[:alpha,:]
    img = plt.imshow(Xapprox, cmap = cm.gray)
    plt.show()