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
# Convert the FF, FF-iPSC+BMP4 or H1+BMP4 MethylC-Seq mapped reads files from the Lister et al. 2011 (Nature) paper (downloaded from http://neomorph.salk.edu/ips_methylomes/data.html on 08/05/2012) to BAM format.
############################################################################################################################################################################################

### TODOs ###
############################################################################################################################################################################################
# TODO: Might add MD and NM tags with samtools calmd
############################################################################################################################################################################################

### INPUT FILE FORMAT ###
############################################################################################################################################################################################
## locations = <integer>; the number of locations the read mapped to (should always be 1).
## mismatches/read-length = NOT USED <integer>; Strictly, the number of mismatches, but equal to the read-length when the quality string is null or absent.
## score = NOT USED <integer>; sum of the quality values of matching bases. Without the quality string this is read length * max quality (ie not useful).
## assembly = <string>; 1, 2, ..., X, Y, M, L. L = lambda spike-in control sample
## strand = <char>; "+" or "-". Which strand the read aligned to.
## position = <integer>; 0-based leftmost Watson-strand mapping coordinate of read.
## name = <string>; read ID.
## copies = <integer>; number of exact copies of the read, i.e. PCR duplicates.
## sequence = <string>; sequence of read. Reverse-complemented if read is aligned to "-" strand.
## quality = <string>; quality scores.
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
#############################################################################################################################################################################################
def makeRNAME(assembly):
    rname = ''.join(['chr', assembly])
    return rname

def makeQNAME(name):
    qname = name
    return qname
    
def makePOS(position):
    start = int(position) # 0-based leftmost mapping position
    return start

def makeFLAG(strand):
    flag = 0x0
    if strand == '-':
        flag += 0x10
    return flag
    
def makeMAPQ():
    return 255 # No mapping information available


def makeCIGAR(sequence):
    cigar = [(0, len(sequence))]
    return cigar

def makeRNEXT():
    return '*'

def makePNEXT():
    return 0

def makeTLEN():
    return 0

def makeSEQ(sequence, strand):
    if strand == '+':
        seq = sequence
    elif strand == '-':
        seq = DNAComplement(sequence[::-1]) # Reverse-complement
    return seq

def DNAComplement(strand):
    return strand.translate(string.maketrans('TAGCNtagcn', 'ATCGNATCGN'))
        
def makeQUAL(quality, sequence, strand):
    if strand == '+':
        qual = quality[:len(sequence)] # If read aligns to OT-strand and quality string is longer than the sequence string, take the "leftmost" quality scores
    elif strand == '-':
        qual = quality[-len(sequence):] # If read aligns to OB-strand and quality string is longer than the sequence string, take the "rightmost" quality scores
    return qual 
    
def makeXG(strand):
    if strand == '+':
        XG = 'CT'
    elif strand == '-':
        XG = 'GA'
    return ('XG', XG)

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
    line = line.rstrip('\n').rsplit('\t')
    # Fields of the Lister-style file
    locations = line[0]
    readlength = line[1]
    score = line[2]
    assembly = line[3]
    strand = line[4]
    position = line[5]
    name = line[6]
    copies = line[7]
    sequence = line[8]
    quality = line[9]
    # Make the SAM/BAM fields
    RNAME = makeRNAME(assembly)
    QNAME = makeQNAME(name)
    FLAG = makeFLAG(strand)
    POS = makePOS(position)
    MAPQ = makeMAPQ()
    CIGAR = makeCIGAR(sequence)
    RNEXT = makeRNEXT() 
    PNEXT= makePNEXT()
    TLEN = makeTLEN()
    SEQ = makeSEQ(sequence, strand)
    QUAL = makeQUAL(quality, sequence, strand)
    XG = makeXG(strand)
    # Single-end 
    read = pysam.AlignedRead()
    read.rname = BAM.gettid(RNAME)
    read.qname = QNAME
    read.flag = FLAG
    read.pos = POS
    read.mapq = MAPQ
    read.cigar = CIGAR
    read.rnext = BAM.gettid(RNEXT)
    read.pnext = PNEXT
    read.tlen = TLEN
    read.seq = SEQ
    read.qual = QUAL
    read.tags = read.tags + [XG]
    if not read.is_paired:
        if read.opt('XG') == 'CT':
            read.tags = read.tags + [('XR', 'CT')]
        elif read.opt('XG') == 'GA':
            read.tags = read.tags + [('XR', 'CT')]
    elif read.is_paired:
        if read.opt('XG') == 'CT' and read.is_read1:
            read.tags = read.tags + [('XR', 'CT')]
        elif read.opt('XG') == 'CT' and read.is_read2:
            read.tags = read.tags + [('XR', 'GA')]
        elif read.opt('XG') == 'GA' and read.is_read1:
            read.tags = read.tags + [('XR', 'CT')]
        elif read.opt('XG') == 'GA' and read.is_read2:
            read.tags = read.tags + [('XR', 'GA')]
    BAM.write(read)
    linecounter += 1
############################################################################################################################################################################################

### Close files
############################################################################################################################################################################################
INFILE.close()
BAM.close()
############################################################################################################################################################################################
