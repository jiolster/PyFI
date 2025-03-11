#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:47:15 2024

@author: joaquin
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import filters, measure, morphology, segmentation, feature
from scipy import ndimage as ndi
from cellpose import models, io
import csv #Export the results as .csv
import errno # for checking whether the montage direcotry already exists



################

# Makes a new direcory for the given path, unless it alredy exists
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise   
            
            
# Makes a montage of the analysed images in each fov
def montage(Nuc, Tc, Tub, segmented_amastigotes, cell_masks):
    
    fig, ax = plt.subplots(nrows=2, ncols=3) #, figsize=(32, 64)
    
    #DAPI
    ax[0][0].imshow(Nuc, cmap='gray')
    ax[0][0].set_title('Hoechst')
    ax[0][0].axis('off')
     
    #PAR o NFkB
    ax[0][1].imshow(Tc, cmap='gray')
    ax[0][1].set_title("anti-T. cruzi")
    ax[0][1].axis('off')
        
    #Tubulina
    ax[0][2].imshow(Tub, cmap='gray')
    ax[0][2].set_title('Tubulina')
    ax[0][2].axis('off')
    
    #Merge
    ax[1][0].imshow(np.moveaxis((Tub, Tc, DAPI), 0, -1))  
    ax[1][0].set_title('Merge')
    ax[1][0].axis('off')
    
    
    #Amastigote masks
    ax[1][1].imshow(segmented_amastigotes)
    ax[1][1].set_title('Amastigotes, %s' %(np.max(segmented_amastigotes)))
    ax[1][1].axis('off')
    
    #Cell masks
    ax[1][2].imshow(cell_masks)
    ax[1][2].set_title('Cells, %s' %(np.max(cell_masks)))
    ax[1][2].axis('off')
    
    plt.tight_layout()
    
    return fig

            
#############


# Wroking direcotry (where program is saved)
wd = "/home/joaquin/Desktop/20250210 - Tul 10 vs 20/Medicion"
os.chdir(wd)

# Montage direcotry
montagedir = os.path.join(wd, "Montage")

make_sure_path_exists(montagedir)


# Folder where all the images are stored, each within a direcotry for all the images in a field
main_folder = '/home/joaquin/Desktop/20250210 - Tul 10 vs 20/Fotos'

# List of each folder containing the images for the fields of view
fields = [f for f in os.listdir(main_folder) if os.path.isdir(os.path.join(main_folder, f))]
    

#Columns for the output table of Average FOV measurments
#List of lists, each row corresponds to the measurments for a field of view
head = ["MOI", "Linea", "Campo", "Celula", "Amastigotes"]

# First row is the name of the columns
results = [head] 


io.logger_setup()
model = models.Cellpose(model_type='cyto3')


for fov in range(len(fields)):
    #List of all images for a field
    campo = os.path.join(main_folder, fields[fov])
    fotos =  [f for f in os.listdir(campo) if os.path.isfile(os.path.join(campo, f))] # Excludes Metadatafolder
    fotos.sort()
    
    # Save the experimental conditions for the fov
    conditions = fields[fov].split("-")
    
    #Load images
    DAPI = plt.imread(os.path.join(campo, fotos[0])) #Nuclear stain
    Tc = plt.imread(os.path.join(campo, fotos[2])) #Trypanosoma cruzi stain
    Tub = plt.imread(os.path.join(campo, fotos[3])) #Cytoplasmatic stain
    
    smooth_amastigotes = filters.gaussian(Tc)
    
    #threshold = filters.threshold_yen(smooth_amastigotes)

    amastigote_mask = smooth_amastigotes > 0.06
    
    amastigote_mask = morphology.remove_small_objects(amastigote_mask, min_size=200)
    
    distance = ndi.distance_transform_edt(amastigote_mask) 
    
    local_max_coords = feature.peak_local_max(distance, min_distance=7, exclude_border=False)
    local_max_mask = np.zeros(distance.shape, dtype=bool)
    local_max_mask[tuple(local_max_coords.T)] = True
    markers = measure.label(local_max_mask)
    
    segmented_amastigotes = segmentation.watershed(-distance, markers, mask=amastigote_mask)
    
    amastigote_props = measure.regionprops(segmented_amastigotes)
    
    
    imgs = [np.moveaxis((Tub, Tc, DAPI), 0, -1)]
    
    # Merge configuration
    # define CHANNELS to run segementation on
    # grayscale=0, R=1, G=2, B=3
    # channels = [cytoplasm, nucleus]
    # if NUCLEUS channel does not exist, set the second channel to 0
    channels = [[0,1]]
    
    if conditions[1] == "Vero":
        diam = 300
        flow = 1
        cellprob = 0
    if conditions[1] == "Caco2":
        diam = 250
        flow = 2
        cellprob = -2
    if conditions[1] == "AC16":
        diam = 300
        flow = 0.5
        cellprob = -1
    if conditions[1] == "HeLa":
         diam = 215
         flow = 0.4
         cellprob = 0
    if conditions[1] == "BeWo":
        diam = 200
        flow = 1
        cellprob = -3
    
    masks, flows, styles, diams = model.eval(imgs, diameter=diam, channels=channels, flow_threshold=flow, cellprob_threshold=cellprob)
    
    clean_masks = morphology.remove_small_objects(masks[0], min_size=5000) #Remove small artifacts
    clean_masks = measure.label(clean_masks) #Rename labels
    cellNum = int(np.max(clean_masks)) #Number of cells detected
    cells = clean_masks.copy() #Image with the masks for all cells
    dim = np.shape(cells) #Dimensions of the image
    
    for c in range(1, cellNum + 1):
        single_cell = cells == c #Choose one of the segmented cells ()
        infection = 0
        for a in range(len(amastigote_props)):
            y, x = amastigote_props[a].centroid
            amastigote_center = (int(y), int(x))
            if single_cell[amastigote_center]:
                infection = infection + 1
    
        #Append measurments to the output table for the whole FOV
        row = [conditions[0], conditions[1], conditions[2], c, infection]
        results.append(row)
    
    # Make montage
    montage_name = "%s.png" % (fields[fov]) 
    montage(DAPI, Tc, Tub, segmented_amastigotes, clean_masks) #Opens a plot with the three images  
    plt.savefig(os.path.join(montagedir, montage_name), bbox_inches='tight', dpi = 300) #Saves the plot
    plt.close() #Closes the plot

#If the output file is not in the program's folder, try at C:/Users/Usuario
with open('results.csv', 'w', newline='') as f: #Measurements for each field's average value
    writer = csv.writer(f)
    writer.writerows(results)
    
'''
vals = Tc.mean(axis=1).flatten()
# plot histogram with 255 bins
b, bins, patches = plt.hist(vals, 255)
plt.xlim([0,255])
plt.show()


fig, ax = filters.try_all_threshold(Tc, figsize=(10, 8), verbose=False)
plt.show()
'''