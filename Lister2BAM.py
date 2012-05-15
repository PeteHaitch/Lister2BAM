#!/usr/bin/python
import pysam
import string

## This program is Copyright (C) 2012, Peter Hickey (hickey@wehi.edu.au)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

### Program description ###
############################################################################################################################################################################################
# Convert the MethylC-Seq mapped reads files from the Lister et al. 2011 (Nature) paper (downloaded from http://neomorph.salk.edu/ips_methylomes/data.html on 08/05/2012-?) to BAM format.
# SequenceA is always interpreted as read1, and sequenceB as read2, for paired-end data.
############################################################################################################################################################################################

### TODOs ###
############################################################################################################################################################################################
# Skip header line
# Check insert size calculations
# Add command line call to @PG tag
# Add line counter
# Don't need the DNA class
# Add arg.parse(), including 'filenames', 'pairedEnd', 'header', 'faidx'
# Make SAM header
# Print error message in make<FIELD>() functions
# Might add MD and NM tags with samtools calmd
# Extend to paire-end data - need to keep track of paired-reads. Keep 'unseen' read-pairs in a dictionary, once seen extract that element (kill its entry) and convert the read-pair. Perhaps read the entire file into a dictionary where keys are read-names, then process the dictionary to create the SAM/BAM. 
############################################################################################################################################################################################

### INPUT FILE FORMAT ###
############################################################################################################################################################################################
## assembly = chromosome name (numeric, hg18)
## strand = strand for which the read is informative ("+" = OT, "-" = OB)
## start = start of read1 (0-based position)
## end =  end of read2 (1-based), i.e. intervals are of the form (start, stop] = {start + 1, start + 2, ..., stop}
## sequenceA = sequence of read1 in left-to-right-Watson-strand orientation. Sequence complemented if strand = "-"
## sequenceB = sequence of read2 (paired-end) or blank (single-end) in left-to-right-Watson-strand orientation. Sequence complemented if strand = "-"
## id = read-ID
# start < end by definition


### Function definitions ###
# All functions return a 2-tuple - the first element for read1 and the second element for read2. If single-end data then only the first element should be used and the second element is set to None
#############################################################################################################################################################################################
def makeRNAME(assembly, sequenceB):
    rname = ''.join(['chr', assembly])
    if sequenceB != '':
        return rname, rname
    else:
        return rname, None

def makeQNAME(RNAME1, ID, sequenceB):
    qname = '_'.join([RNAME1, ID])
    if sequenceB != '':
        return qname, qname
    else:
        return qname, None
    
def makePOS(start, end, sequenceA, sequenceB):
    if sequenceB != '': # Read is paired-end
        start1 = int(start) # 0-based leftmost mapping position
        start2 = int(end) - len(sequenceB) # 0-based leftmost mapping position
        return start1, start2
    else:
        start = int(start) # 0-based leftmost mapping position
        return start, None
        
def makeFLAG(sequenceA, sequenceB, strand):
    if sequenceB != '': # Read is paired-end in-sequencing
        flag1 = 0x0 # Flag value for read1 in a paired-end read
        flag2 = 0x0 # Flag value for read2 in a paired-end read
        if sequenceB != '': # Read is paired-end in-sequencing
            flag1 += 0x01
            flag2 += 0x01
            flag1 += 0x02 # Flag is properly-paired according to the aligner (forcing to be true)
            flag2 += 0x02 # Flag is properly-paired according to the aligner (forcing to be true)
        if strand == '+':
            flag1 += 0x20
            flag2 += 0x10
        elif strand == '-':
            flag1 += 0x10
            flag2 += 0x20
        else:
            print 'No strand set for read', RNAME1            
        flag1 += 0x40 # Assuming sequenceA refers to read1 in the read-pair
        flag2 += 0x80 # Assuming sequenceB refers to read2 in the read-pair
        return flag1, flag2
    else: # Read is single-end
        flag = 0x0
        if strand == '-':
            flag += 0x10
        return flag, None
    
def makeMAPQ(sequenceB):
    if sequenceB != '':
        return 255, 255 # No mapping information available
    else:
        return 255, None # No mapping information available

def makeCIGAR(sequenceA, sequenceB):
    if sequenceB != '':
        cigar1 = [(0, len(sequenceA))]
        cigar2 = [(0, len(sequenceB))]
        return cigar1, cigar2
    else:
        cigar = [(0, len(sequenceA))]
        return cigar, None

def makeRNEXT(RNAME1, RNAME2):
    if RNAME2 is not None:
        return RNAME2, RNAME1
    else:
        return '*', None

def makePNEXT(start1, start2):
    if start2 is not None:
        return start2, start1
    else:
        return 0, None

def makeTLEN(start, end, sequenceB):
    if sequenceB != '': # Paired-end read
        abs_tlen = int(end) - int(start) # absolute value of TLEN 
        return abs_tlen, -abs_tlen
    else:
        return 0, None

def makeSEQ(sequenceA, sequenceB, strand):
    if strand == '+':
        seq1 = sequenceA
        seq2 = sequenceB 
    elif strand == '-':
        seq1 = sequenceA
        seq2 = sequenceB
    if sequenceB != '':
        return seq1, seq2
    else:
        return seq1, None
        
def makeQUAL(sequenceA, sequenceB):
    qual1 = 'E' * len(sequenceA)
    qual2 = 'E' * len(sequenceB)
    if sequenceB != '':
        return qual1, qual2
    else:
        return qual1, None

def createHeader():
    FAIDX = open('/export/share/disk501/lab0605/Bioinformatics/databases/bahlolab_db/WGBS/genomes/hg18+lambda_phage/hg18+lambda.fa.fai', 'r')
    faidx = FAIDX.read().rstrip().rsplit('\n')
    hd = {'VN': '1.0', 'SO': 'unsorted'}
    sq = []
    for i in range(0, len(faidx)):
        line = faidx[i].rsplit('\t')
        sq.append({'LN': int(line[1]), 'SN': line[0]})
    pgid = 'Lister2BAM.py'
    vn = '1.0'
    cl = 'command line'
    pg = [{'ID': pgid, 'VN': vn, 'CL': cl}]
    header = {'HD': hd, 'SQ': sq, 'PG': pg}
    FAIDX.close()
    return header


#############################################################################################################################################################################################
# Make SAM header
header = createHeader()

# Make BAM file
BAM = pysam.Samfile('test.bam', 'wh', header = header) # This writes SAM rather than BAM

# Loop over files chromosome-by-chromosome
FILENAMES = ['../Lister_2011_BS-seq_data/ADS/methylCseq_reads_ads/methylCseq_reads_ads_10']

for FILENAME in FILENAMES:
    f = open(FILENAME, 'r')
    linecounter = 1
    for line in f: # Loop over the gzip file line-by-line and convert to corresponding SAM/BAM entry
        if linecounter == 1: # Skip the header line
            linecounter +=1
            continue
        line = line.rsplit()
        # Fields of the Lister-style file
        assembly = line[0]
        strand = line[1]
        start = line[2]
        end = line[3]
        sequenceA = line[4]
        sequenceB = line[5]
        ID = line[6]
        # Make the SAM/BAM fields
        RNAME1, RNAME2 = makeRNAME(assembly, sequenceB)
        QNAME1, QNAME2 = makeQNAME(RNAME1, ID, sequenceB)
        FLAG1, FLAG2 = makeFLAG(sequenceA, sequenceB, strand)
        POS1, POS2 = makePOS(start, end, sequenceA, sequenceB)
        MAPQ1, MAPQ2 = makeMAPQ(sequenceB)
        CIGAR1, CIGAR2 = makeCIGAR(sequenceA, sequenceB)
        RNEXT1, RNEXT2 = makeRNEXT(RNAME1, RNAME2) 
        PNEXT1, PNEXT2 = makePNEXT(POS1, POS2)
        TLEN1, TLEN2 = makeTLEN(start, end, sequenceB)
        SEQ1, SEQ2 = makeSEQ(sequenceA, sequenceB, strand)
        QUAL1, QUAL2 = makeQUAL(sequenceA, sequenceB)
        if sequenceB != '': # Paired-end read
            read1 = pysam.AlignedRead()
            read2 = pysam.AlignedRead()
            read1.rname = BAM.gettid(RNAME1)
            read2.rname = BAM.gettid(RNAME2)
            read1.qname = QNAME1
            read2.qname = QNAME2
            read1.flag = FLAG1
            read2.flag = FLAG2
            read1.pos = POS1
            read2.pos = POS2
            read1.mapq = MAPQ1
            read2.mapq = MAPQ2
            read1.cigar = CIGAR1
            read2.cigar = CIGAR2
            read1.rnext = BAM.gettid(RNEXT1)
            read2.rnext = BAM.gettid(RNEXT2)
            read1.pnext = PNEXT1
            read2.pnext = PNEXT2
            read1.tlen = TLEN1
            read2.tlen = TLEN2
            read1.seq = SEQ1
            read2.seq = SEQ2
            read1.qual = QUAL1
            read2.qual = QUAL2
            BAM.write(read1)
            BAM.write(read2)
        elif sequenceB == '': # Single-end 
            read = pysam.AlignedRead()
            read.rname = BAM.gettid(RNAME1)
            read.qname = QNAME1
            read.flag = FLAG1
            read.pos = POS1
            read.mapq = MAPQ1
            read.cigar = CIGAR1
            read.rnext = BAM.gettid(RNEXT1)
            read.pnext = PNEXT1
            read.tlen = TLEN1
            read.seq = SEQ1
            read.qual = QUAL1
            OUT.write(read)
        linecounter += 1
    f.close()

BAM.close()
