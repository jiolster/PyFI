# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 16:32:16 2024

@author: Joaquin
"""
# Libraries


#Set system path to image

import os

folder = 'F:/13-7-23 - BeWo Activacion de PAR Tul wt/Fotos'
image_name = 'Tul_BeWo_0h_1_1_PAR.tif'
filepath = os.path.join(folder, image_name)

#Read image

from matplotlib import image as mpimg

image = mpimg.imread(filepath) #image as ndarray

#Convert image to grayscale

import numpy as np

# FUNCTION 1
# Dot product between image as np.array and RGB weights
# Returns a 64bit grayscale image
def rgb2gray(rgb):
    grayscale = np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])
    return grayscale


#Otsu threshold
from skimage.filters import threshold_otsu

# FUNCTION 2
# Makes a binary mask (True or False for each pixel) for a grayscale image using Otsu's method.
# If the whole cell is stained, the function marks as True those pixels that form the cells.
# DAPI staining sets nuclei to True.

def otsu_mask(image):
    thresh = threshold_otsu(image) #Value 
    binary = image > thresh
    return binary


# FUNCTION 3
# Set background (threshold >) to NaN, only retains the fluorescence emitted by the cells 
# Requires a grayscale image from the fluorescent channel that needs to be measured and 
# a binary mask separating cells from background
# Returns the segmented image as an array

def background2nan(image, mask):
    cell = image.copy()
    dim = np.shape(cell)
    for i in range(dim[0]):
        for j in range(dim[1]):
            if not mask[i, j]:
                cell[i, j] = np.nan
    return cell


# FUNCTION 4
# Requires an image for whole cells with the background set to NaN, and a binary mask for de nuclei
# Returns two images: the nuclear fluorescence and the citosolic fluorescence

def segment_cell(cells, nuclei):
    nuclear = cells.copy()
    cytosolic = cells.copy()
    dim = np.shape(cells)
    for i in range(dim[0]):
        for j in range(dim[1]):
            if nuclei[i, j]:
                cytosolic[i, j] = np.nan
            if not nuclei[i, j]:
                nuclear[i, j] = np.nan     
    return nuclear, cytosolic


###Auto load images
    
from os.path import isfile, join 

files =  [f for f in os.listdir(folder) if isfile(join(folder, f))] # Makes a list of the names of all the images in the folder
'''
#Makes a list of all the conditions included in an image´s name (separated by "_")
conditions = files[0].split('_') #Makes a list of all the conditions included in an image´s name (separated by "_")
parasite_strain = conditions[0] # What strain was used in the experiment
cell_line = conditions[1] # What cell line was infected during the experiment
Treatment = conditions[2] # How long the cells were infected, or what else was done to them
Sample = conditions[3] # The number for the physical sample that was taken with the above conditions
Field = conditions[4] # The sample's field from which the image was obtained
Channel = conditions[5] # The stain that was observed when aquiring the image
'''
files.sort() # Sort the image names so that images from each field are adjacent each other within the list
    
from itertools import groupby

# Makes a list of lists. Each sublist is made up of the conditions for each image.
field_identifiers = sorted([file.split("_") for file in files])


# List of lists. Each sublist of lists with the keywords for each image in a field
grouped = [list(value) for key, value in groupby(field_identifiers, lambda x: x[:-1])] 

# List of lists. All images for each field are grouped in a sublist. 
# Re-Joins the keywords for each image, recomposing the original names. 
fields = []
for group in grouped:
    temp = []
    for i in range(len(group)):
        temp.append("_".join(group[i]))
    fields.append(temp)


export = [["Cepa", "Celula", "Tratamiento", "Muestra", "Campo", "PAR total", "PAR nuclear", "PAR citoplastamico", "Relacion"]]


for i in range(len(fields)):
    filepathDAPI = os.path.join(folder, fields[i][1])
    DAPI = mpimg.imread(filepathDAPI)
    DAPI = rgb2gray(DAPI) 

    filepathPAR = os.path.join(folder, fields[i][2])
    PAR = mpimg.imread(filepathPAR)
    PAR = rgb2gray(PAR) 

    filepathTc = os.path.join(folder, fields[i][3])
    Tc = mpimg.imread(filepathTc)
    Tc = rgb2gray(Tc) 
    
    # Mask for segmenting cells   
    contorno = otsu_mask(Tc)

    # Channel of interest with the background set to NaN
    signal = background2nan(PAR, contorno)

    # Mask for segmenting nuclei
    nucelos = otsu_mask(DAPI)

    # Channel of interest segmented into nuclear and citoplasmatic signal
    nuc, citosol = segment_cell(signal, nucelos)

    # Overall Mean signal intensity
    total = np.nanmean(signal)

    #Nuclear mean signal intensity
    nuclear = np.nanmean(nuc)

    #Citoplasmatic signal intensity
    citoplasmatico = np.nanmean(citosol)
    
    # Nuclear / citoplasmatic signal
    nc = nuclear/citoplasmatico
    
    row = [grouped[i][0][0], grouped[i][0][1], grouped[i][0][2], grouped[i][0][3], grouped[i][0][4], total, nuclear, citoplasmatico, nc]
    
    export.append(row)
    
    

import csv

#Saving at C:/Users/Usuario
with open('out.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(export)

# Convert array to PIL Image
from PIL import Image

img = Image.fromarray(total)

# Save as TIFF image
img.save("output.tif")

