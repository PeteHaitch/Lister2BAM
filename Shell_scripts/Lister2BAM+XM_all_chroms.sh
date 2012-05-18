# Peter Hickey
# 18/05/2012
# A shell script to convert all aligned reads from a MethylC-seq experiment from Lister 2009 or 2011 from Lister-style to BAM format.
# There is generally a separate Lister-style file for each chromosome so we process each chromosome seaprately in pseudo-parallel, and then merge the results
# Script should be run in folder containing the Lister-style files, i.e. methylC_seq_reads_<sample_ID>_<chrom>

# Sample-specific variables
SAMPLE_ID=
LISTER_PREFIX=

# Make directory for BAM files
mkdir ../BAM/

# Fixed variables (shouldn't need changing)
LISTER2BAM='python /home/users/lab0605/hickey/Lister2BAM/Lister2BAM.py'
REF=/export/share/disk501/lab0605/Bioinformatics/databases/bahlolab_db/WGBS/genomes/hg18+lambda_phage/hg18+lambda.fa
REF_INDEX=${REF}.fai
XM_TAG='python /home/users/lab0605/hickey/Lister2BAM/XM_tag.py'
BAM_DIR=../BAM
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M L)

for i in ${chroms[@]}
do
FILE=${LISTER_PREFIX}${i}
${LISTER2BAM} ${FILE} ${BAM_DIR}/${SAMPLE_ID}_${i}.bam ${REF_INDEX} &
done

wait

for i in ${chroms[@]}
do
FILE=${LISTER_PREFIX}${i}
${XM_TAG} ${BAM_DIR}/${SAMPLE_ID}_${i}.bam ${BAM_DIR}/XM_${SAMPLE_ID}_${i}.bam ${REF} &
done

wait