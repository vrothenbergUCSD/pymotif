To test PyMotif on experimental data, please run the following command (adjust output directory and number of cores as needed):
````
pymotif experimental-data/H3K4me3_peaks.txt experimental-data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa ./yeastOutput/ -size 200 -mask -cores 9
````

Optionally, you can reduce runtime by modifying the number of peaks to analyze:
````
pymotif experimental-data/H3K4me3_peaks.txt experimental-data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa ./yeastOutput/ -size 200 -mask -cores 9 -num_peaks 2000
````

To compare results to HOMER, please run the following command (after properly installing and configuring HOMER):
````
findMotifsGenome.pl experimental-data/H3K4me3_peaks.txt experimental-data/Saccharomyces_cerevisiae.R64-1-1.forHomer.fa ./yeastOutput_homer/ -mask -size 200
 ````
 
 ### Sources
 `H3K4me3_peaks.txt` is from Gene Expression Omnibus (GEO) at NCBI (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2590466). You can manually download it with
 ````
 wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2590nnn/GSM2590466/suppl/GSM2590466_1H3K4me3_WT_1_peaks.txt.gz
 ````
 `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa` is from Ensembl (https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/). You can manually download it with 
 ````
 wget https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
 gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
 ````
  `Saccharomyces_cerevisiae.R64-1-1.forHomer.fa` is the same file as above, but with modified headings for each chromosome as required for running HOMER. Specifically, headings were modified to include "chr" before the chromosome number. You can manually recreate this file with
   ````
 wget https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
 gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
 sed -i 's/>/>chr/g' Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
 ````
