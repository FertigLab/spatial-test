import os
import subprocess as sp
import shutil
import h5py
import csv
import pyarrow as pa
import pyarrow.parquet as pq

# Take paramaters from params.csv
param_file = open('params.csv','r',newline='')
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

if hd:
    input_path = input_path + '/binned_outputs/' + resolution

# List of all files that the pipeline uses
files = ['filtered_feature_bc_matrix.h5',
         'spatial/tissue_positions_list.csv',
         'spatial/scalefactors_json.json',
         'spatial/tissue_hires_image.png',
         'spatial/tissue_lowres_image.png']
if hd:
    files[1] = 'spatial/tissue_positions.parquet'
    files.append('../../feature_slice.h5')
    files.append('../../molecule_info.h5')

# Ensure the needed files are present
for file in files:
    assert os.path.exists(input_path+'/'+file)

# Create output folders
if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.mkdir(output_path)
if hd:
    os.mkdir(output_path + '/binned_outputs')
    os.mkdir(output_path + '/binned_outputs/'+resolution)
    os.mkdir(output_path + '/binned_outputs/'+resolution+'/spatial')
    output_path = output_path + '/binned_outputs/' + resolution
else:
    os.mkdir(output_path + '/spatial')
    # If not HD, prepare to read the original csv
    old_spots = open(input_path+'/'+files[1],'r',newline='')
    reader = csv.reader(old_spots, delimiter=',')

# Copy files to output folders
for file in files:
    shutil.copyfile(input_path+'/'+file,output_path+'/'+file)
os.chdir(output_path)

# Open tissue position list and filtered feature matrix
filt_mat = h5py.File(files[0],'r+')
if hd:
    positions = pq.read_table(files[1])
else:
    os.remove(files[1])
    new_spots = open(files[1],'w',newline='')
    writer = csv.writer(new_spots, delimiter=',')

# Filter spots from the position list
codes = []
if hd:
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
    print('Parquet file processed')
else:
    for row in reader:
        if row[1] == '1':
            if xmin<=int(row[3])<=xmax and ymin<=int(row[2])<=ymax:
                writer.writerow(row)
                codes.append(row[0])
    print("CSV file processed")
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
filt_mat['matrix'].create_dataset('data',data = counts)
filt_mat['matrix'].create_dataset('indices',data = inds)
filt_mat['matrix/shape'][1] = len(spots)

print('Size reduced to '+str(len(spots))+' spots')
if len(spots)<20:
    print("Warning: area is too small for pipeline to process")
elif len(spots)<200:
    print('Warning: area may be too small for Spacemarkers to run properly')

filt_mat.close()
if not hd:
    old_spots.close()
    new_spots.close()

# Repack the h5 file to save space
shutil.copyfile(files[0],'temp_'+files[0])
os.remove(files[0])
sp.run(['h5repack','temp_'+files[0],files[0]])
os.remove('temp_'+files[0])