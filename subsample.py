import os
import sys
import subprocess as sp
import shutil
import h5py
import csv
import json
import pyarrow as pa
import pyarrow.parquet as pq

# Take paramaters from params.csv
if len(sys.argv) != 2:
    print("Usage: python subsample.py <params.csv>")
    sys.exit(1)

params_arg = sys.argv[1]

param_file = open(params_arg,'r',newline='')
param_read = csv.reader(param_file, delimiter=',')
assert param_read.__next__() == ['input_path','output_path','hd_resolution','xmin','xmax','ymin','ymax']
params = param_read.__next__()
param_file.close()

input_path = params[0]
output_path = params[1]
resolution = params[2]
hd = bool(len(resolution))
xmin = int(params[3])
xmax = int(params[4])
ymin = int(params[5])
ymax = int(params[6])


# List of all files that the pipeline uses
files = ['filtered_feature_bc_matrix.h5',
         'spatial/tissue_positions_list.csv',
         'spatial/scalefactors_json.json',
         'spatial/tissue_hires_image.png',
         'spatial/tissue_lowres_image.png']
if hd:
    files[1] = 'spatial/tissue_positions.parquet'
    base_path = input_path
    input_path = input_path + '/binned_outputs/' + resolution

# Ensure the needed files are present
for file in files:
    print(f"Checking for {file} in {input_path}")
    assert os.path.exists(input_path+'/'+file)

# Create output folders and create dummy feature_slice.h5
if os.path.exists(output_path):
    shutil.rmtree(output_path)

if hd:
    os.makedirs(output_path + '/binned_outputs/'+resolution+'/spatial')
    slice_src = os.path.join(base_path, 'feature_slice.h5')
    slice_dst = os.path.join(output_path, 'feature_slice.h5')
    shutil.copyfile(slice_src, slice_dst)

    #slim down feature_slice.h5 to keep size minimal
    with h5py.File(slice_dst, 'r+') as f:
        objects = list(f.keys())
        for obj in objects:
            del f[obj]
    sp.run(['h5repack',slice_dst,os.path.join(output_path,'temp_feature_slice.h5')])
    shutil.move(os.path.join(output_path,'temp_feature_slice.h5'),slice_dst)
    output_path = output_path + '/binned_outputs/' + resolution
else:
    os.makedirs(output_path + '/spatial')


# Copy files to output folders
for file in files:
    src = os.path.join(input_path, file)
    dst = os.path.join(output_path, file)
    print(f'Copying {src} to {dst}')
    shutil.copyfile(src,dst)
os.chdir(output_path)


# Filter spots from the position list
codes = []
if hd:
    positions = pq.read_table(files[1])
    print(f"Original dataset: {pa.compute.max(positions[3])} by {pa.compute.max(positions[2])} bins")
    spot_table = []
    for i in range(positions.num_columns):
        spot_table.append([])

    for i in range(positions.num_rows):
        if positions[1][i].as_py() == 1:
            if xmin<=positions[3][i].as_py()<=xmax and ymin<=positions[2][i].as_py()<=ymax:
                for j in range(positions.num_columns):
                    spot_table[j].append(positions[j][i].as_py())
                codes.append(positions[0][i].as_py())
    table = pa.table({'barcode': spot_table[0],
                    'in_tissue': spot_table[1],
                    'array_row': spot_table[2],
                    'array_col': spot_table[3],
                    'pxl_row_in_fullres': spot_table[4],
                    'pxl_col_in_fullres': spot_table[5]})
    os.remove(files[1])
    pq.write_table(table, files[1])
    print('Parquet file processed')
else:
    os.remove(files[1])
    new_spots = open(files[1],'w',newline='')
    writer = csv.writer(new_spots, delimiter=',')
    old_spots = open(os.path.join(input_path,files[1]),'r',newline='')
    reader = csv.reader(old_spots, delimiter=',')
    
    xvals = []
    yvals = []
    for row in reader:
        if row[1] == '1':                   # Only keep spots in tissue
            xvals.append(int(row[3]))
            yvals.append(int(row[2]))
            if xmin<=int(row[3])<=xmax and ymin<=int(row[2])<=ymax:
                writer.writerow(row)
                codes.append(row[0])
    print(f"Original dataset: {max(xvals)} by {max(yvals)} spots")
    print(f"CSV file processed")
    old_spots.close()
    new_spots.close()

codes.sort()

# Create reduced feature matrix
filt_mat = h5py.File(files[0],'r+')
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
filt_mat['matrix'].create_dataset('data',data = counts, compression='gzip')
filt_mat['matrix'].create_dataset('indices',data = inds, compression='gzip')
filt_mat['matrix/shape'][1] = len(spots)
filt_mat.close()

print('Size reduced to '+str(len(spots))+' spots')
if len(spots)<20:
    print("Warning: area is too small for pipeline to process")
elif len(spots)<200:
    print('Warning: area may be too small for Spacemarkers to run properly')


# Repack the h5 file to save space
sp.run(['h5repack',files[0],'temp_'+files[0]])
shutil.move('temp_'+files[0], files[0])

# Pretend that lowres is the highres image and update scalefactors_json.json
shutil.copyfile('spatial/tissue_lowres_image.png', 'spatial/tissue_hires_image.png')
# Update scalefactors_json.json
with open('spatial/scalefactors_json.json', 'r') as f:
    scalefactors = json.load(f)
scalefactors['tissue_hires_scalef'] = scalefactors['tissue_lowres_scalef']
with open('spatial/scalefactors_json.json', 'w') as f:
    json.dump(scalefactors, f, indent=4)

