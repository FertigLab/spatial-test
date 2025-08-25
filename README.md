The subsample program takes a Visium or Visium HD dataset and reduces its size to improve the pipeline's runtime for testing purposes. It takes an input dataset and makes a copy containing only the files necessary to the pipeline, then reduces the size of the files containing the data.

The program takes its parameters from params.csv, which must be in the same directory as subsample.py.
There are seven parameters: input_path, output_path, hd_resolution, xmin, xmax, ymin, ymax

input_path is the directory path of the original dataset that will have its size reduced. The program does not alter any files in this directory. This path can be relative, such as "../../sample_data", or explicit, such as "/Users/iuy/.nextflow/assets/break-through-cancer/btc-spatial-pipelines/sample_data". For explicit paths, start the path with a "/" symbol.

output_path is the path where the new, smaller dataset will be created. This path be relative or explicit in the same way as input_path.
--------
WARNING: If files or directories already exist in this location, they will be deleted.
--------

hd_resolution is the resolution the program will act on when working on an HD dataset. Examples: "square_002um", "square_008um", "square_016um". For regular Visium datasets, leave this value empty.

The last four parameters define the rectangular area which is kept to make the smaller subsample. Note that some datasets have irregular dimensions. The program will print the dimensions of the original input dataset for reference.