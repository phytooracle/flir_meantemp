# Flir mean temperature extractor

Extracts the mean temperature of plants from plots.
 
## Inputs

Directory of plotclipped flir files

## Outputs

Single csv file containing the mean temperature of plots and plants within a plot

## Arguments and Flags
- **Required Arguments:** 
    - **Directory of plotclipped flirs:** 'dir'
    - **Date of scan date:** '-d', 'scan_date'
    - **geojson shapefile:** '-g'

- **Optional Arguments**
    - **Output directory:** '-o', '--outdir', default='meantemp_out/'
                                        
## Adapting the Script

The script requires a geojson file of the field. Such file needs to contain plot information. The script also needs to modied in line `138`: the input directory is **two** directories above the tif files; modify to suit your organizational structure.


##Executing example (using singularity)
`singularity run -B $(pwd):/mnt --pwd /mnt/ docker://phytooracle/flir_meantemp -g <shp.geojson> -o <out_dir/> -d <scan_date> <tif_dir>`
