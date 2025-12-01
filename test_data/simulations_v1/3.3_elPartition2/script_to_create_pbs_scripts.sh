#!/bin/bash

# Base file to copy
template="3.2_t10_len100.pbs"

# Loop over all t and len combinations
for t in 10 20 40 80 160; do
  for len in 100 1000 10000; do
    newfile="3.3_t${t}_len${len}.pbs"
    newscript="3.3_t${t}_len${len}_script.pg"

    # Copy the template to the new file
    cp "$template" "$newfile"

    # Replace the script name inside the file
    sed -i "s/3.2_t10_len100_script.pg/${newscript}/g" "$newfile"
  done
done

