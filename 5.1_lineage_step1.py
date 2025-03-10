## The mitochondrial mutation sites of different cells were obtained by comparing mitochondrial reference sequences.
# Code operation mode： python ../puxi_step1.py ./ ./C230517002-chr_q30.bam ./PBMC2_Plasma_cellname.txt PBMC2 Plasma

import pysam
import os
import sys

# Mitochondrial sequences were extracted from bam files of single cells
os.system("samtools view -h anno_decon.bam chrM >anno_decon_chrM.bam") 
# The sequences matching q30 were screened
os.system("samtools view -@ 10 -h -q 30 anno_decon_chrM.bam -o anno_decon_chrM_q30.bam")

mydir = sys.argv[1]  # Working directory
rawbam = sys.argv[2]  # bam file name  ../anno_decon_chrM_q30.bam
selectcell = sys.argv[3] # txt file name, which contains only a list of single-cell names, eg：PVAT8_B_cellname.txt
prename = sys.argv[4]  # Output file name prefix, eg:"PVAT8"
celltype = sys.argv[5]  # Cell type name， eg: "B"

os.chdir(mydir)
os.getcwd()

## Reads the name of the specified single cell
bc_set = set()
with open(selectcell) as barcodes:
    for (index, line) in enumerate(barcodes):
        bc = line.strip()
        bc_set.add(bc)

bam = pysam.AlignmentFile(rawbam,"r")


## Screen for sequences corresponding to cell names
pairedreads = pysam.AlignmentFile("{0}_{1}_chrM_q30.bam".format(prename,celltype), "wb", template=bam)
list1 = list()
for (index,read) in enumerate(bam):
    if read.get_tag("DB") in bc_set:
        pairedreads.write(read)
    if index % 3000000 == 0:
        print("have been processed",index)

print("end")
pairedreads.close()


# Sort by "DB" name
os.system("samtools sort -t DB {0}_{1}_chrM_q30.bam  -o {0}_{1}_chrM_q30_sorted.bam".format(prename,celltype))


### Python 3.6.8
## Gets a single bam file for the selected cell
import pysam

### Input varibles to set
# file to split on
unsplit_file = "{0}_{1}_chrM_q30_sorted.bam".format(prename,celltype)
# where to place output files
out_dir = "./{0}_{1}_chrM_bam_split/".format(prename,celltype)
os.mkdir(out_dir)

# variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile(unsplit_file, "rb")
for read in samfile.fetch( until_eof=True):
    # barcode itr for current read
    CB_itr = read.get_tag('DB')
    # if change in barcode or first line; open new file  
    if( CB_itr!=CB_hold or itr==0):
        # close previous split file, only if not first read in file
        if( itr!=0):
            split_file.close()
        CB_hold = CB_itr
        itr+=1
        split_file = pysam.AlignmentFile( out_dir + "{}_{}.bam".format(prename,CB_hold), "wb", template=samfile)
        # print(itr)

    # write read with same barcode to file
    split_file.write(read)
print("split cell end")
split_file.close()
samfile.close()



os.chdir(out_dir)
os.getcwd()

## Build an index for each single-cell bam file
os.system("for i in `ls *.bam`; do samtools index $i; done")  

## Get the ATCG variation data for each cell bam file, stored in the "processed.MAE_mito.rds" file, and ppl2_run.py can get from https://github.com/songjiajia2018/ppl
os.system("python ./ppl2_run.py -t 10 --input-filelist -p -m -r --name singlecell --input PVAT8_B_cellname.txt --outdir ../")