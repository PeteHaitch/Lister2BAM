#!/usr/bin/python
import pysam
import string
import argparse

# The MIT License (MIT)

# Copyright (c) [2014] [Peter Hickey (peter.hickey@gmail.com)]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

### Program description ###
############################################################################################################################################################################################
# Convert the ADS, ADS-adipose or ADS-iPSC MethylC-Seq mapped reads files from the Lister et al. 2011 (Nature) paper (downloaded from http://neomorph.salk.edu/ips_methylomes/data.html on 08/05/2012) to BAM format.
############################################################################################################################################################################################

### TODOs ###
############################################################################################################################################################################################
# TODO: Might add MD and NM tags with samtools calmd
############################################################################################################################################################################################

### INPUT FILE FORMAT ###
############################################################################################################################################################################################
## assembly = chromosome name (numeric, hg18)
## strand = strand for which the read is informative ("+" = OT, "-" = OB)
## start = start of read1 (0-based position)
## end =  end of read2 (1-based), i.e. intervals are of the form (start, stop] = {start + 1, start + 2, ..., stop}
## sequenceA = sequence of read1 in left-to-right-Watson-strand orientation. Sequence complemented if strand = "-"
## sequenceB = sequence of read2 in left-to-right-Watson-strand orientation. Sequence complemented if strand = "-"
## id = read-ID
# NB: start < end by definition
############################################################################################################################################################################################

### Command line passer ###
############################################################################################################################################################################################
parser = argparse.ArgumentParser(description='Convert Lister-style alignment files of MethylC-Seq data to BAM format.')
parser.add_argument('infile', metavar = 'infile',
                   help='The filename of the Lister-style file that is to be converted to BAM format')
parser.add_argument('outfile', metavar = 'out.bam',
                  help='The path to the new SAM/BAM file.')
parser.add_argument('ref_index', metavar = 'reference.fa.fai',
                  help='The path to the index (.fai file) of reference genome FASTA file.')
args = parser.parse_args()
############################################################################################################################################################################################

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
    return strand.translate(string.maketrans('TAGCNtagcn', 'ATCGNATCGN'))
        
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
        return ('XG', XG), None

def createHeader():
    FAIDX = open(args.ref_index, 'r')
    faidx = FAIDX.read().rstrip().rsplit('\n')
    hd = {'VN': '1.0', 'SO': 'unsorted'}
    sq = []
    for i in range(0, len(faidx)):
        line = faidx[i].rsplit('\t')
        sq.append({'LN': int(line[1]), 'SN': line[0], 'AS': 'hg18+lambda_phage'})
    pgid = 'Lister_style_6_to_bam.py'
    vn = '1.0'
    cl = ' '.join([pgid, args.infile, args.outfile, args.ref_index])
    pg = [{'ID': pgid, 'VN': vn, 'CL': cl}]
    header = {'HD': hd, 'SQ': sq, 'PG': pg}
    FAIDX.close()
    return header
#############################################################################################################################################################################################

### Open files ###
############################################################################################################################################################################################
INFILE = open(args.infile, 'r')
header = createHeader()
BAM = pysam.Samfile(args.outfile, 'wb', header = header)
############################################################################################################################################################################################

### The main loop ###
############################################################################################################################################################################################
# Loop over methylC_seq_reads files file-by-file (i.e. chromosome-by-chromosome)
print 'Input file is', args.infile
linecounter = 1
for line in INFILE: # Loop over the file line-by-line and convert to an AlignedRead instance
    if linecounter == 1: # Skip the header line
        linecounter +=1
        continue
    line = line.rstrip('\n').rsplit('\t')
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
    if sequenceA == '':
        print 'WARNING: Empty sequenceA at line', linecounter, 'in file', args.infile
        print line
    if sequenceB == '':
        print 'WARNING: Empty sequenceB at line', linecounter, 'in file', args.infile
        print line
    # Paired-end: using readL/readR notation, thus for the Lister protocol a OT-strand readL=read1 and readR=read2 whereas for OB-strand readL=read2 and readR=read1
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
    if not readL.is_paired:
        if readL.opt('XG') == 'CT':
            readL.tags = readL.tags + [('XR', 'CT')]
        elif readL.opt('XG') == 'GA':
            readL.tags = readL.tags + [('XR', 'CT')]
    elif readL.is_paired:
        if readL.opt('XG') == 'CT' and readL.is_readL1:
            readL.tags = readL.tags + [('XR', 'CT')]
        elif readL.opt('XG') == 'CT' and readL.is_readL2:
            readL.tags = readL.tags + [('XR', 'GA')]
        elif readL.opt('XG') == 'GA' and readL.is_readL1:
            readL.tags = readL.tags + [('XR', 'CT')]
        elif readL.opt('XG') == 'GA' and readL.is_readL2:
            readL.tags = readL.tags + [('XR', 'GA')]
    if not readR.is_paired:
        if readR.opt('XG') == 'CT':
            readR.tags = readR.tags + [('XR', 'CT')]
        elif readR.opt('XG') == 'GA':
            readR.tags = readR.tags + [('XR', 'CT')]
    elif readR.is_paired:
        if readR.opt('XG') == 'CT' and readR.is_readR1:
            readR.tags = readR.tags + [('XR', 'CT')]
        elif readR.opt('XG') == 'CT' and readR.is_readR2:
            readR.tags = readR.tags + [('XR', 'GA')]
        elif readR.opt('XG') == 'GA' and readR.is_readR1:
            readR.tags = readR.tags + [('XR', 'CT')]
        elif readR.opt('XG') == 'GA' and readR.is_readR2:
            readR.tags = readR.tags + [('XR', 'GA')]
    BAM.write(readL)
    BAM.write(readR)
    linecounter += 1
############################################################################################################################################################################################

### Close files
############################################################################################################################################################################################
INFILE.close()
BAM.close()
############################################################################################################################################################################################
