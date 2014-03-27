Python scripts to convert Lister-style alignment files from Lister et al. Nature (2009 and 2011) to `BAM` format (Binary Sequence Alignment/Map format, see http://samtools.sourceforge.net/SAM1.pdf). 

The Lister-style files can be downloaded from http://neomorph.salk.edu/human_methylome/data.html (2009 data) and http://neomorph.salk.edu/ips_methylomes/ (2011 data).

__WARNING__: While I have attempted to make the output files compatible with Bismark's `SAM` format so that, for example, they y nice with `bismark_methylation_extractor`, there are a couple of problems. Please see the issue at https://github.com/PeteHaitch/Lister2BAM/issues/1 for a patch to `bismark_methylation_extractor` which fixes two problems.
