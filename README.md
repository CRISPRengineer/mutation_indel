# mutation_indel
This script will analyze sam alignment files to calculate the total insertions and deletions (indels) for a set of files within a given directory. 

./mutation.pl -i <indir> -b <left bound> -l length
  
  -i: the input directory, default cwd.
  -b: the left boundary to begin searching for indels
  -l: the length of the window to search for indels, default = 30.
  
  The output will be displayed in the terminal and written to a file called summary.txt. Here is an example data summary:

sample	total reads	mapped reads	deletion reads	insertion reads	PE indel reads	PE indel freq
ABEgR1_ABE_hNeg3	2025	2025	2	0	2	0.0987654320987654
ABEgR1_BE4-gam_ABEgR1	3802	3802	72	1	72	1.89374013677012
ABEgR1_Cas9_ABEgR1	3761	3761	1869	405	2260	60.0904014889657
ABEgR1_ABE_ABEgR1	50	50	0	0	0	0

*PE indels are only mutations that are confirmed by both R1 and R2 for highest accuracy and sensitivity of rare mutations.
