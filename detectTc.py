#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:47:15 2024

@author: joaquin
"""

from skimage.feature import shape_index


from skimage.filters import threshold_yen
def yen_mask(image, blur = False, erode = 0, fill = 0, remove_small = 0):
    if blur:
        image = gaussian(image)
    thresh = threshold_yen(image) #Value 
    binary = image > thresh
    if erode != 0:
        binary = morphology.isotropic_erosion(binary, erode)
    if fill != 0:
        binary = morphology.isotropic_closing(binary, fill)
    if remove_small > 0:
        binary = morphology.remove_small_objects(binary,  min_size = remove_small)
    return binary




kinetos = yen_mask(DAPI)
plt.imshow(kinetos)
HighTub = yen_mask(Tub, erode = 1, fill = 5, remove_small=100)
plt.imshow(HighTub)
mitosis = yen_mask(Tub, erode = 1, fill = 5, remove_small=500)
plt.imshow(mitosis)
amas = HighTub ^ mitosis
plt.imshow(amas)
numAmas = countCells(amas)
print(numAmas)
amas1 = label(amas)
plt.imshow(amas1)
amas2 = regionprops(amas1)
print(int(amas2[0].centroid[0]))
np.max(amas1)

plt.imshow(kinetos * amas1)
