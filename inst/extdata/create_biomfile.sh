#!/bin/bash

# this script converts the soil_euk tsv reads table into a biom file
#
# first, install the biom-format python package (http://biom-format.org/index.html):
# pip install numpy
# pip install biom-format
# pip install h5py

# convert tsv into biom
biom convert -i data-raw/litiere_euk_reads_4biom.tsv -o data-raw/litiere_euk_reads_hdf5.biom --table-type="OTU table" --to-hdf5
# delete tsv file
rm data-raw/litiere_euk_reads_4biom.tsv