import os
import shutil
import numpy as np
import h5py
import csv

# Parameters that should probably be external
user = 'iuy'
old_data = 'sample_data'
new_data = 'sample_small'
xmin = 150
xmax = 250
ymin = 150
ymax = 250

# File paths for the old and new datasets
parent = '/Users/' + user + '/.nextflow/assets/break-through-cancer/btc-spatial-pipelines/'
sample_path = parent + old_data + '/'
new_path = parent + new_data + '/'

# List of all files that the pipeline uses
files = ['filtered_feature_bc_matrix.h5',
         'spatial/tissue_positions_list.csv',
         'spatial/scalefactors_json.json',
         'spatial/tissue_hires_image.png',
         'spatial/tissue_lowres_image.png']

if os.path.exists(parent + 'out_small'):
    shutil.rmtree(parent + 'out_small')

# Open tissue position list
os.chdir(sample_path)
for file in files:
    assert os.path.exists(file)
old_spots = open(files[1],'r',newline='')
reader = csv.reader(old_spots, delimiter=',')

# Copy old dataset to new path
if os.path.exists(new_path):
    shutil.rmtree(new_path)
os.mkdir(new_path)
os.mkdir(new_path + 'spatial')
for file in files:
    shutil.copyfile(sample_path+file,new_path+file)

# Open new tissue position list and filtered feature matrix
os.chdir(new_path)
os.remove(files[1])
new_spots = open(files[1],'w',newline='')
writer = csv.writer(new_spots, delimiter=',')
filt_mat = h5py.File(files[0],'r+')

# Copy only certain rows from the csv
codes = []

for row in reader:
    if row[1] == '1':
        if xmin<=int(row[3])<=xmax and ymin<=int(row[2])<=ymax:
            writer.writerow(row)
            codes.append(row[0].encode())
codes.sort()

# Create lists of reduced data size
spots = []
ptrs = [0]
counts = []
inds = []
total = 0
barcodes = filt_mat['matrix/barcodes'][:].tolist()
indptr = filt_mat['matrix/indptr'][:]
indices = filt_mat['matrix/indices'][:]
datas = filt_mat['matrix/data'][:]

for code in codes:
    spot = barcodes.index(code)
    imin = indptr[spot]
    imax = indptr[spot+1]
    count = datas[imin:imax]
    ind = indices[imin:imax]

    spots.append(spot)
    for i in range(imax-imin):
        counts.append(count[i])
        inds.append(ind[i])
        total += 1
    ptrs.append(total)

# Delete h5 datasets and replace with reduced lists
del filt_mat['matrix/barcodes']
del filt_mat['matrix/data']
del filt_mat['matrix/indices']
del filt_mat['matrix/indptr']
filt_mat['matrix'].create_dataset('barcodes',data = codes)
filt_mat['matrix'].create_dataset('indptr',data = ptrs)
filt_mat['matrix'].create_dataset('data',data = counts)
filt_mat['matrix'].create_dataset('indices',data = inds)
filt_mat['matrix/shape'][1] = len(spots)

print('Size reduced to '+str(len(spots))+' spots')
old_spots.close()
new_spots.close()
filt_mat.close()