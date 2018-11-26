# excise column of sample run ids for downloading
cut -d "," -f 1 data/metadata/SraRunTable.txt | tail -n +2 | \
    head -n -1 > data/metadata/accessions_to_download.txt

# download target sra files from the file list using aspera ascp
# documentation here:
# https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=prefetch
prefetch.2.8.2 -t \
  ascp -a "/usr/local/bin/ascp|/mnt/raid/class_data/2017-03-Fall/downloading_li2015/misc/asperaweb_id_dsa.openssh" \
  -p 0.5 \
  --option-file data/metadata/accessions_to_download.txt
