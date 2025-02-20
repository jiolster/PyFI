#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:43:27 2024

@author: joaquin
"""

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
def montage(Nuc, Tc, Tub, segmented_amastigotes):
    
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
    
    plt.tight_layout()
    
    return fig

            
#############


# Wroking direcotry (where program is saved)
wd = "/home/joaquin/Desktop/MOI 1/Medicion"
os.chdir(wd)

# Montage direcotry
montagedir = os.path.join(wd, "Montage")

make_sure_path_exists(montagedir)


# Folder where all the images are stored, each within a direcotry for all the images in a field
main_folder = '/home/joaquin/Desktop/MOI 1/Fotos'

# List of each folder containing the images for the fields of view
fields = [f for f in os.listdir(main_folder) if os.path.isdir(os.path.join(main_folder, f))]
    

#Columns for the output table of Average FOV measurments
#List of lists, each row corresponds to the measurments for a field of view
head = ["Cepa", "Celula", "Tiempo", "Muestra", "Campo", "Threshold", "Objetos"]

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
    
    smooth_amastigotes = filters.gaussian(Tc)
    
    threshold = filters.threshold_otsu(smooth_amastigotes)

    amastigote_mask = smooth_amastigotes > threshold
    
    amastigote_mask = morphology.remove_small_objects(amastigote_mask, min_size=200)
     
    distance = ndi.distance_transform_edt(amastigote_mask) 
    
    local_max_coords = feature.peak_local_max(distance, min_distance=7, exclude_border=False)
    local_max_mask = np.zeros(distance.shape, dtype=bool)
    local_max_mask[tuple(local_max_coords.T)] = True
    markers = measure.label(local_max_mask)
    
    segmented_amastigotes = segmentation.watershed(-distance, markers, mask=amastigote_mask)
    
    
    #Append measurments to the output table for the whole FOV
    row = [conditions[0], conditions[1], conditions[2], conditions[3], conditions[4], threshold, np.max(segmented_amastigotes)]
    results.append(row)
    
    
    # Make montage
    montage_name = "%s.png" % (fields[fov]) 
    montage(DAPI, Tc, Tub, segmented_amastigotes) #Opens a plot with the three images  
    plt.savefig(os.path.join(montagedir, montage_name), bbox_inches='tight', dpi = 300) #Saves the plot
    plt.close() #Closes the plot
    
#If the output file is not in the program's folder, try at C:/Users/Usuario
with open('results.csv', 'w', newline='') as f: #Measurements for each field's average value
    writer = csv.writer(f)
    writer.writerows(results)
  
x = []
y = []
for i in range(2, len(results)):
    x.append(results[i][5])
    y.append(results[i][6])
    

fig, ax = plt.subplots()
ax.scatter(x, y, s = 3)

ax.set(ylim=(0, 150), xticks=np.arange(0, 0.4, 0.03))

plt.show()


vals = Tc.mean(axis=1).flatten()
# plot histogram with 255 bins
b, bins, patches = plt.hist(vals, 255)
plt.xlim([0,50])
plt.show()


fig, ax = filters.try_all_threshold(Tc, figsize=(10, 8), verbose=False)
plt.show()
