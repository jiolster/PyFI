# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:31:58 2024

@author: Joaquin
"""

#Libraries

import os #Read directories
from os.path import isfile, isdir, join #Identify files, make paths
from matplotlib import image as mpimg #Load images as arrays
import matplotlib.pyplot as plt # For making the montages
import cmasher as cmr # For better colourmaps
import numpy as np #Manipulate arrays
from skimage import filters #Thresholding and guassian filter
from skimage.filters import threshold_otsu
from skimage.filters import threshold_minimum
from skimage.filters import gaussian
from skimage.filters import laplace
from skimage.measure import label, regionprops #For counting
import csv #Export the results as .csv
import errno # for checking whether the montage direcotry already exists

###Functions

# FUNCTION 1
# Dot product between image as np.array and RGB weights
# Returns a 64bit grayscale image
def rgb2gray(rgb):
    grayscale = np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])
    return grayscale


# FUNCTION 2
# From a direcotry with all the images from a field, 
# find the slice with the best focus by nuclear brightness
# Change channel to focus with by changing "fc = ch##.tif"
# Leica software exports images using different indeces if the number of slices is over or under 10
# specify number of channels to get number of slices
def brightest(fotos, channels, fc = "ch00.tif"):
    b = []
    for i in range(len(fotos)):
        foto = os.path.join(campo, fotos[i])
        args = fotos[i].split("_")
        if args[-1] == fc:
            image = mpimg.imread(foto)
            brightness = np.sum(image)
            b.append(brightness)   
            
    if (len(fotos) / channels) > 10:
        if b.index(max(b)) > 9:
            z = "z" + str(b.index(max(b))) #Brightest slice
        else:
            z = "z0" + str(b.index(max(b))) #Brightest slice
    else:
        z = "z" + str(b.index(max(b))) #Brightest slice
        
    return z   
   
# FUNCTION 3
# Find best focus by nuclear sharpness
# Change channel to focus with by changing "fc = ch##.tif"
# Leica software exports images using different indeces if the number of slices is over or under 10
# specify number of channels to get number of slices
def sharpest(fotos, channels, fc = "ch00.tif"):
    s = []
    for i in range(len(fotos)):
        foto = os.path.join(campo, fotos[i])
        args = fotos[i].split("_")
        if args[-1] == fc:
            image = mpimg.imread(foto)
            gaus = gaussian(image)
            filtered = laplace(gaus)
            sharpness = np.sum(filtered)
            s.append(sharpness)
            
    if (len(fotos) / channels) > 10:
        if s.index(max(s)) > 9:
            z = "z" + str(s.index(max(s))) #Sharpest slice
        else:
            z = "z0" + str(s.index(max(s))) #Sharpest slice
    else:
        z = "z" + str(s.index(max(s))) #Sharpest slice
        
    return z

# FUNCTION 4
# List of the images for each channel in a field, at the focal plane (fov).
# If the fov isn't a z stack, returns the images for each channel.
# If you want to exclude bright field, specify its channel number BF = "ch##.tif".

def inFocus(fotos, focal, fov, BF = "None"):
    focus = []
    if BF != "None":
        for i in range(len(fotos)):
            args = fotos[i].split("_")
            if len(args) > 2:
                if args[-2] == focal and args[-1] != BF:
                    foto = os.path.join(campo, fotos[i])
                    focus.append(foto)
            else:
                if args[-1] != BF:
                    foto = os.path.join(campo, fotos[i])
                    focus.append(foto)
    else:
        for i in range(len(fotos)):
            args = fotos[i].split("_")
            if len (args) > 2:
                if args[-2] == focal:
                    foto = os.path.join(campo, fotos[i])
                    focus.append(foto)
            else:
                foto = os.path.join(campo, fotos[i])
                focus.append(foto)
    return focus

# FUNCTION 5
# Makes a binary mask (True or False for each pixel) for a grayscale image using Otsu's method.
# If the whole cell is stained, the function marks as True those pixels that form the cells.
# DAPI staining sets nuclei to True.

def otsu_mask(image, blur = False):
    if blur:
        image = gaussian(image)
    thresh = threshold_otsu(image) #Value 
    binary = image > thresh
    return binary

# FUNCTION 6
# Binary mask using triangle method. Does not remove objects on borders
def triangle_mask(image, blur = False):
    if blur:
        image = gaussian(image) # Takes care of holes, but fills in dim images 
    thresh = filters.threshold_triangle(image) # Value 
    binary = image > thresh #All pixels above the value return True forming a binary mask
    return binary

# FUNCTION 7
# Binary mask using minimum method. Does not remove objects on borders

def minimum_mask(image, blur = False):
    if blur:
        image = gaussian(image) # Takes care of holes, but fills in dim images 
    thresh = filters.threshold_minimum(image) # Value 
    binary = image > thresh #All pixels above the value return True forming a binary mask
    return binary


# FUNCTION 8
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

# FUNCTION 9
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

# FUNCTION 10
# Counts number of objects in a binary mask.
def countCells(mask):
    count = 0
    label_image = label(mask)
    for region in regionprops(label_image):
        # take regions with large enough areas
        if region.area >= 200:
            count += 1
    return count
    

#FUNCTION 11
# Makes a new direcory for the given path, unless it alredy exists
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise   
            
#Function 12
# Makes a montage of the analysed images in each fov
def montage(Nuc, Measure, Tub, field):
    
    image_names = field.split("-")
    
    fig, ax = plt.subplots(ncols=3, figsize=(32, 32))
    
    #DAPI
    ax[0].imshow(DAPI, cmap='cmr.jungle')
    ax[0].set_title('DAPI', fontsize = 36)
    ax[0].axis('off')
     
    #PAR o NFkB
    ax[1].imshow(Measure, cmap='cmr.rainforest')
    ax[1].set_title(image_names[-1], fontsize = 36)
    ax[1].axis('off')
        
    #Tubulina
    ax[2].imshow(Tub, cmap='afmhot')
    ax[2].set_title('Tubulina', fontsize = 36)
    ax[2].axis('off')
    
    return fig



################

# Wroking direcotry (where program is saved)
wd = "C:\\Users\\gaizq\\Desktop\\20240809 - Vero Tul PAR y NFkB\\Medicion"
os.chdir(wd)

# Montage direcotry
montagedir = join(wd, "Montage")

make_sure_path_exists(montagedir)


# Folder where all the images are stored, each within a direcotry for all the images in a field
main_folder = 'C:\\Users\\gaizq\\Desktop\\20240809 - Vero Tul PAR y NFkB\\Fotos'

# List of each folder containing the images for the fields of view
fields = [f for f in os.listdir(main_folder) if isdir(join(main_folder, f))]



#Columns for the output table
head = ["Tratamiento", "Celula", "Condicion", "Muestra", "Campo", "Tincion", "Total", "Nuclear", "Citoplasmatico", "Rel","Max", "Nucleos"]

#List of lists, each row corresponds to the measurments for a field of view

# First row is the name of the columns
results = [head] 

# Loop over the direcotries for each field of view 
for i in range(len(fields)):
    #List of all images for a field
    campo = os.path.join(main_folder, fields[i])
    fotos =  [f for f in os.listdir(campo) if isfile(os.path.join(campo, f))] # Excludes Metadatafolder
    
    # Select the images for the most in focus plane
    z = brightest(fotos, 3)
    focus = inFocus(fotos, z, campo) # 0 = DAPI, 2 = PAR, 1 = Tub
    
    #Load images
    DAPI = mpimg.imread(focus[0]) #Nuclear stain
    Measure = mpimg.imread(focus[2]) #Image of interest
    Tub = mpimg.imread(focus[1]) #Cytoplasmatic stain
    
    #Masks 
    nuc = otsu_mask(DAPI) #Nuclei
    cp = triangle_mask(Tub, blur = True) #Cytoplasm
    
    # Set background to NaN
    Measure = np.float_(Measure) #Convert 8 bit to 64bit float
    MeasureNaN = background2nan(Measure, cp)
    
    #Separate different parts of the image
    nuclear, cytoplasmatic = segment_cell(MeasureNaN, nuc)
    
    #Count nuclei in image
    countNuc = countCells(nuc)
    
    #Measure mean fluorescence intensity
    totalMean = np.nanmean(MeasureNaN)
    nucMean = np.nanmean(nuclear)
    cpMean = np.nanmean(cytoplasmatic)
    rel = nucMean / cpMean 
    
    norm = np.max(Measure) #Save the maximum pixel intensity for normalization
    
    # Save the field's name
    conditions = fields[i].split("-")
    
    #Append measurments to output table.
    row = [conditions[0], conditions[1], conditions[2], conditions[3], conditions[4], conditions[5], totalMean, nucMean, cpMean, rel, norm, countNuc]
    results.append(row)
    
    
    # Make montage
    montage_name = "%s.png" % (fields[i]) 
    montage(DAPI, Measure, Tub, fields[i]) #Opens a plot with the three images  
    plt.savefig(join(montagedir, montage_name), bbox_inches='tight') #Saves the plot
    plt.close() #Closes the plot

#If the output file is not in the program's folder, try at C:/Users/Usuario
with open('results.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(results)
    