import os
import shutil
import numpy as np
import h5py
import pyarrow as pa
import pyarrow.parquet as pq

# Parameters that should probably be external
user = 'iuy'
sample_folder = 'sample_hd_2'
resolution = 'square_016um'
xmin = 71
xmax = 190
ymin = 0
ymax = 85
genes = 18085

# File paths for the dataset
parent = '/Users/' + user + '/.nextflow/assets/break-through-cancer/btc-spatial-pipelines/'
sample_path = parent + sample_folder + '/binned_outputs/' + resolution + '/'

# List of all files that the pipeline uses
files = ['filtered_feature_bc_matrix.h5',
         'spatial/tissue_positions.parquet',
         'spatial/scalefactors_json.json',
         'spatial/tissue_hires_image.png',
         'spatial/tissue_lowres_image.png']

# Ensure the needed files are present
os.chdir(sample_path)
for file in files:
    assert os.path.exists(file)

# Create backups of the original files
backups = ['filtered_feature_bc_matrix_full.h5',
           'spatial/tissue_positions_full.parquet']
for i in [0,1]:
    if not os.path.exists(sample_path + backups[i]):
        shutil.copyfile(sample_path+files[i],sample_path+backups[i])
    else:
        os.remove(sample_path+files[i])
        shutil.copyfile(sample_path+backups[i],sample_path+files[i])

# Open tissue position list and filtered feature matrix
filt_mat = h5py.File(files[0],'r+')
positions = pq.read_table(files[1])

# Copy only certain rows from the parquet
codes = []
spot_table = []
for i in range(positions.num_columns):
    spot_table.append([])

for i in range(positions.num_rows):
    if positions[1][i].as_py() == 1:
        if xmin<=int(positions[3][i])<=xmax and ymin<=int(positions[2][i])<=ymax:
            for j in range(positions.num_columns):
                spot_table[j].append(positions[j][i])
            codes.append(positions[0][i].as_py())
table = pa.table({'barcode': spot_table[0],
                  'in_tissue': spot_table[1],
                  'array_row': spot_table[2],
                  'array_col': spot_table[3],
                  'pxl_row_in_fullres': spot_table[4],
                  'pxl_col_in_fullres': spot_table[5]})
os.remove(files[1])
pq.write_table(table, files[1])
codes.sort()

print('Parquet file processed')

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
    spot = barcodes.index(code.encode())
    imin = indptr[spot]
    imax = indptr[spot+1]
    count = datas[imin:imax]
    ind = indices[imin:imax]

    spots.append(spot)
    for i in range(imax-imin):
        if ind[i] < genes:
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
filt_mat['matrix/shape'][0] = genes
filt_mat['matrix/shape'][1] = len(spots)

filt_mat['matrix/features/feature_type'].resize((genes,))
filt_mat['matrix/features/genome'].resize((genes,))
filt_mat['matrix/features/id'].resize((genes,))
filt_mat['matrix/features/name'].resize((genes,))

print('Size reduced to '+str(len(spots))+' spots and '+str(genes)+' genes')
filt_mat.close()