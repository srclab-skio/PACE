#!/bin/bash

# Specify the directory containing the NetCDF files
input_directory="/Volumes/SRC_HDD_Mas/Data/PACE_Sapelo_Matchups"
output_directory="//Volumes/SRC_HDD_Mas/Data/PACE_Sapelo_Matchups/PACE_Sapelo_Matchups"

# Define the polygon for subsetting (same for all files)
geoRegion='POLYGON ((-82.17597 27.275438, -82.17597 33.2286567, -75.834220886 33.2286567, -75.02487182 33.2286567, -82.17597 27.275438))'

# Loop over all NetCDF files in the directory
for input_file in "$input_directory"/*.nc; do
  # Get the base filename without the directory path
  base_filename=$(basename "$input_file" .nc)

  # Define the output filename by appending '_subset' to the original filename
  output_file="${output_directory}/${base_filename}_subset.nc"

  # Run the GPT Subset command for each file
  /Applications/snap/bin/gpt Subset -Ssource="$input_file" -PgeoRegion="$geoRegion" -t "$output_file" -f NetCDF4-CF

  # Print a message after processing each file
  echo "Subsetted $input_file and saved to $output_file"
done