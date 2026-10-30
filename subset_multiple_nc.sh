#!/bin/bash

# Specify the directory containing the NetCDF files
input_directory="/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_3_2/Irrawady"
output_directory="/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_3_2/Irrawady"

# Define the polygon for subsetting (same for all files)
 # geoRegion='POLYGON ((-82.17597 27.275438, -82.17597 33.2286567, -75.834220886 33.2286567, -75.02487182 33.2286567, -82.17597 27.275438))' # SAB
  # geoRegion='POLYGON ((-82.17597 27.275438, -82.17597 34, -75.834220886 34, -75.02487182 34, -82.17597 27.275438))' # New SAB
# geoRegion='POLYGON ((-81.6 30.5, -81.6 32, -80.5 32, -80.5 30.5, -81.6 30.5))' # Altamaha River
# geoRegion='POLYGON ((83 18.4,87.5 22.5 , 92 23, 94.7 18.4, 83 18.4))' # BoB
# geoRegion='POLYGON ((86 19,86 23 , 92.5 23, 92.5 19, 86 19))' # BoB rivers
 # geoRegion='POLYGON ((116 37, 121 42 , 127 35, 120.7 32.35, 116 37))' # Yellow sea
# geoRegion='POLYGON ((118 37, 118 39 , 120 39, 120 37, 118 37))' # Yellow river
# geoRegion='POLYGON ((-51.6 4, -46.2 4.5, -46 -2, -51.4 -2,-51.6 4))' # Amazon River
 geoRegion='POLYGON ((94 14.5,94 16.3 , 96 16.3, 96. 14.5, 94 14.5))' # Irrawaddy River rivers

# Loop over all NetCDF files in the directory
for input_file in "$input_directory"/*.nc; do
  # Get the base filename without the directory path
  base_filename=$(basename "$input_file" .nc)

  # Define the output filename by appending '_subset' to the original filename
  output_file="${output_directory}/${base_filename}_subset.nc"

  # Run the GPT Subset command for each file
  /Applications/esa-snap/bin/gpt Subset -Ssource="$input_file" -PgeoRegion="$geoRegion" -t "$output_file" -f NetCDF4-CF

  # Print a message after processing each file
  echo "Subsetted $input_file and saved to $output_file"
done