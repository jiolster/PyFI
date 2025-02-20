#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 17:00:08 2025

@author: joaquin
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import filters, measure, morphology, segmentation, feature
from scipy import ndimage as ndi
import csv #Export the results as .csv
import errno # for checking whether the montage direcotry already exists
import math


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
    
    
    #smooth_DAPI = filters.gaussian(DAPI, sigma = 1.5)
    
    #fig, ax = try_all_threshold(smooth_DAPI, figsize=(10, 8), verbose=False)
    #plt.show()
    
    DAPI_threshold = filters.threshold_li(DAPI)
    DAPI_mask = DAPI > DAPI_threshold
    
    DAPI_mask = ndi.binary_fill_holes(DAPI_mask)
    
    DAPI_mask = morphology.erosion(DAPI_mask)
    DAPI_mask = morphology.remove_small_objects(DAPI_mask, min_size=50)
    
    
    DAPI_mask = segmentation.clear_border(DAPI_mask)
    
    #Watershed
    distance = ndi.distance_transform_edt(DAPI_mask) 
    
    local_max_coords = feature.peak_local_max(distance, min_distance=7, exclude_border=False)
    local_max_mask = np.zeros(distance.shape, dtype=bool)
    local_max_mask[tuple(local_max_coords.T)] = True
    markers = measure.label(local_max_mask)
    
    DAPI_watershed = segmentation.watershed(-distance, markers, mask=DAPI_mask)
    
    DAPI_clean = segmentation.clear_border(DAPI_watershed) 
    
    Nuclei_mask = morphology.remove_small_objects(DAPI_clean, min_size=500)
    
    Kinetoplast_mask = DAPI_clean ^ Nuclei_mask
    
    Nuclei_labels = measure.label(Nuclei_mask)  
    
    Kinetoplast_labels = measure.label(Kinetoplast_mask)
    
    Kinetoplast_props = measure.regionprops(Kinetoplast_labels)
    
    Nuclei_props = measure.regionprops(Nuclei_labels)
    
    #Measure distances between Nuceli centroids and kinetoplast centroids
    
    infected_cells = []
    for tc in range(len(Kinetoplast_props)):
        Kinetoplast_coord = Kinetoplast_props[tc].centroid #(y, x)
    
        yk = int(Kinetoplast_coord[0]) #Kinetoplast y coordinate
        xk = int(Kinetoplast_coord[1]) #Kinetoplast x coordinate
    
    
        distances = []
        for i in range(len(Nuclei_props)):
            Nuclei_coord = Nuclei_props[i].centroid #(y, x)
            yn = int(Nuclei_coord[0]) #Nucleus y coordinate
            xn = int(Nuclei_coord[1]) #Nucleus x coordinate
        
        
            d = math.sqrt((xn - xk)**2 + (yn - yk)**2)
            
            distances.append(d)
    
        infected_cells.append(distances.index(min(distances)))
    
    for cell in range(1, len(Nuclei_props) + 1):
        amastigotes = 0
        for j in range(len(infected_cells)):
            if infected_cells[j] == cell:
                amastigotes = amastigotes + 1
        row = [conditions[0], conditions[1], conditions[2], cell, amastigotes]
        results.append(row)
                
    
    # Make montage
    montage_name = "%s.png" % (fields[fov]) 
    montage(DAPI, Tc, Tub, Kinetoplast_labels, Nuclei_labels) #Opens a plot with the three images  
    plt.savefig(os.path.join(montagedir, montage_name), bbox_inches='tight', dpi = 300) #Saves the plot
    plt.close() #Closes the plot

#If the output file is not in the program's folder, try at C:/Users/Usuario
with open('results.csv', 'w', newline='') as f: #Measurements for each field's average value
    writer = csv.writer(f)
    writer.writerows(results)
        
