## Create a genome browser from GenBank file
import genomebrowser as GB
import urllib.request
# import D3GB, urllib.request

# Download GenBank file
gbk, headers = urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/297/395/GCF_000297395.2_ASM29739v2/GCF_000297395.2_ASM29739v2_genomic.gbff.gz')
# gbk, headers = urllib.request.urlretrieve('file:///home/wing/workspace/gb/Sequence_Files/uex_short.gb')

# Genome browser generation.
# It creates a genome browser ready to view in a Firefox browser.
# For a server version, ready to be shared with Apache as a Website, set server=True
gb = GB.gbk2genomebrowser(gbk,directory='Micromonospora_Lupini_gbk')