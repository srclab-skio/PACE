#!/bin/bash

# Specify the filename prefix
prefix=pixexl

# Loop through the files and generate .txt files, echo is like print
# the filenames will be different for each of the sensors and level. So, you will have to change the naming accordingly
for i in *L2.OC.nc; do 
  echo "Creating file: ${prefix}${i}.txt"
  # Run PixEx with appropriate arguments
original_filename="$i"
  echo $original_filename
  # string_to_remove="SEAHAWK1_HAWKEYE"
  string_to_remove="" # don't need to remove anything for this case

  echo $string_to_remove 
  new_filename=$(echo "$original_filename" | sed "s/$string_to_remove//")
  echo $new_filename
  # SH = 130m = 7x7, MA = 1km = 1 x 1, VIIRS = 500m = 2 x 2, but 3x3 as a requirement, uneven number
  string1="1by1"
  string2="P33M"
  string3="3by3_Mean"
  string4="7by7_Mean"
 
  string5="P55"
  /Applications/snap/bin/gpt PixEx -PcoordinatesFile=./Pins_AG.txt -PsourceProductPaths="$i" -PoutputDir=./Data_extractions/ -PoutputFilePrefix="${string1}${original_filename}"
  # /Applications/snap/bin/gpt; change this to your own path where gpt file exists 
  # Pins_AG.txt	is the pins or lat lon coordinates
  # if you consider only single pixel extraction, you can comment the following lines for mean.
  /Applications/snap/bin/gpt PixEx -PcoordinatesFile=./Pins_AG.txt -PsourceProductPaths="$i" -PoutputDir=./Data_extractions/ -PwindowSize=3 -PoutputFilePrefix="${string3}${original_filename}" -PaggregatorStrategyType='mean'
  /Applications/snap/bin/gpt PixEx -PcoordinatesFile=./Pins_AG.txt -PsourceProductPaths="$i" -PoutputDir=./Data_extractions/ -PwindowSize=7 -PoutputFilePrefix="${string4}${original_filename}" -PaggregatorStrategyType='mean'

  # You can also uncomment the following line if you want to create content files
  # echo "File content for ${prefix}${i}" > "${prefix}${i}.txt"
done

