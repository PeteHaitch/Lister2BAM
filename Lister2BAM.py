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
# Fix all functions so stranded-ness is properly handled
# Add XG tag to reads
# All reads informative for the OB strand have the incorrect orientation when viewed in IGV. 
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
# All functions return a 2-tuple - the first element for readL (the leftmost read, regardless of strand) and the second element for readR (the rightmost read, regardless of strand)
# If single-end data then only the first element should be used and the second element is set to None
#############################################################################################################################################################################################
def makeRNAME(assembly, sequenceB):
    rname = ''.join(['chr', assembly])
    if sequenceB != '':
        rnameL = rname
        rnameR = rname
        return rnameL, rnameR
    else:
        return rname, None

def makeQNAME(RNAMEL, ID, sequenceB):
    qname = '_'.join([RNAMEL, ID])
    if sequenceB != '':
        qnameL = qname
        qnameR = qname
        return qnameL, qnameR
    else:
        return qname, None
    
def makePOS(start, end, sequenceA, sequenceB):
    if sequenceB != '': # Read is paired-end
        startL = int(start) # 0-based leftmost mapping position
        startR = int(end) - len(sequenceB) # 0-based leftmost mapping position
        return startL, startR
    else:
        start = int(start) # 0-based leftmost mapping position
        return start, None
        
def makeFLAG(sequenceA, sequenceB, strand):
    if sequenceB != '': # Read is paired-end in-sequencing
        flagL = 0x0 # Flag value for read1 in a paired-end read
        flagR = 0x0 # Flag value for read2 in a paired-end read
        if sequenceB != '': # Read is paired-end in-sequencing
            flagL += 0x01 # Paired-end read
            flagR += 0x01 # Paired-end read
            flagL += 0x02 # Flag is properly-paired according to the aligner (forcing to be true)
            flagR += 0x02 # Flag is properly-paired according to the aligner (forcing to be true)
            if strand == '+':
                flagL += 0x20 # Seq of readR is reverse-complemented
                flagR += 0x10 # Seq of readR is reverse-complemented
                flagL += 0x40 # Leftmost read is read1
                flagR += 0x80 # Rightmost read is read2
            elif strand == '-':
                flagL += 0x20 # Seq of read1 is reverse-complemented
                flagR += 0x10 # Seq of read1 is reverse-complemented
                flagR += 0x40 # Rightmost read is read1
                flagL += 0x80 # Leftmost read is read2
        else:
            print 'No strand set for read', RNAME1
        return flagL, flagR
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
        cigarL = [(0, len(sequenceA))]
        cigarR = [(0, len(sequenceB))]
        return cigarL, cigarR
    else:
        cigar = [(0, len(sequenceA))]
        return cigar, None

def makeRNEXT(RNAMEL, RNAMER):
    if RNAMER is not None:
        return RNAMER, RNAMEL
    else:
        return '*', None

def makePNEXT(startL, startR):
    if startR is not None:
        return startR, startL
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
        seqL = sequenceA
        seqR = sequenceB 
    elif strand == '-':
        seqL = DNAComplement(sequenceA)
        seqR = DNAComplement(sequenceB)
    if sequenceB != '':
        return seqL, seqR
    else:
        return seqL, None

def DNAComplement(strand):
    return strand.translate(string.maketrans('TAGCtagc', 'ATCGATCG'))
        
def makeQUAL(sequenceA, sequenceB):
    qualL = 'E' * len(sequenceA)
    qualR = 'E' * len(sequenceB)
    if sequenceB != '':
        return qualL, qualR
    else:
        return qualL, None
    
def makeXG(sequenceB, strand):
    if strand == '+':
        XG = 'CT'
    elif strand == '-':
        XG = 'GA'
    if sequenceB != '':
        XGL = ('XG', XG)
        XGR = ('XG', XG)
        return XGL, XGR
    else:
        return ('XG', XG)

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
        RNAMEL, RNAMER = makeRNAME(assembly, sequenceB)
        QNAMEL, QNAMER = makeQNAME(RNAMEL, ID, sequenceB)
        FLAGL, FLAGR = makeFLAG(sequenceA, sequenceB, strand)
        POSL, POSR = makePOS(start, end, sequenceA, sequenceB)
        MAPQL, MAPQR = makeMAPQ(sequenceB)
        CIGARL, CIGARR = makeCIGAR(sequenceA, sequenceB)
        RNEXTL, RNEXTR = makeRNEXT(RNAMEL, RNAMER) 
        PNEXTL, PNEXTR = makePNEXT(POSL, POSR)
        TLENL, TLENR = makeTLEN(start, end, sequenceB)
        SEQL, SEQR = makeSEQ(sequenceA, sequenceB, strand)
        QUALL, QUALR = makeQUAL(sequenceA, sequenceB)
        XGL, XGR = makeXG(sequenceB, strand)
        if sequenceB != '': # Paired-end read - using readL/readR notation, thus for the Lister protocol a OT-strand readL = read1 and readR = read2 whereas OB-strand readL = read2 and readR = read1
            readL = pysam.AlignedRead()
            readR = pysam.AlignedRead()
            readL.rname = BAM.gettid(RNAMEL)
            readR.rname = BAM.gettid(RNAMER)
            readL.qname = QNAMEL
            readR.qname = QNAMER
            readL.flag = FLAGL
            readR.flag = FLAGR
            readL.pos = POSL
            readR.pos = POSR
            readL.mapq = MAPQL
            readR.mapq = MAPQR
            readL.cigar = CIGARL
            readR.cigar = CIGARR
            readL.rnext = BAM.gettid(RNEXTL)
            readR.rnext = BAM.gettid(RNEXTR)
            readL.pnext = PNEXTL
            readR.pnext = PNEXTR
            readL.tlen = TLENL
            readR.tlen = TLENR
            readL.seq = SEQL
            readR.seq = SEQR
            readL.qual = QUALL
            readR.qual = QUALR
            readL.tags = readL.tags + [XGL]
            readR.tags = readR.tags + [XGL]
            BAM.write(readL)
            BAM.write(readR)
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
            read.tags = read.tags + [XGL]
            OUT.write(read)
        linecounter += 1
    f.close()

BAM.close()
