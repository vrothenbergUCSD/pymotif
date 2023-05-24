# Design

The design of PyMotif takes heavy inspiration from the motif analysis tool [HOMER](http://homer.ucsd.edu/homer/motif/) developed primarily by Chris Benner 

## Preprocessing

1. Extraction of Sequences
   1. Appropriate genomic DNA is extracted from background genomic regions e.g. hg19
2. Background Selection
   1. If the background sequences were not explicitly defined, PyMotif will automatically select them.  If you are using genomic positions, sequences will be randomly selected from the genome, matched for GC% content (to make GC normalization easier in the next step).  Custom backgrounds can be specified with "-bg <file>". 
3. GC Normalization
   1. Optional - Sequences in the target and background sets are then binned based on their GC-content (5% intervals).  Background sequences are weighted to resemble the same GC-content distribution observed in the target sequences.  This helps avoid simply finding motifs that are GC-rich when analyzing sequences from CpG Islands. 
4. Autonormalization 
   1. Optional - Often the target sequences have an imbalance in the sequence content other than GC%.  This can be caused by biological phenomenon, such as codon-bias in exons, or experimental bias caused by preferrential sequencing of A-rich stretches etc.  If these sources of bias are strong enough, HOMER will lock on to them as features that significanly differentiate the target and background sequences.  HOMER now offers autonormalization as a technique to remove (or partially remove) imblances in short oligo sequences (i.e. AA) by assigning weights to background sequences.  The proceedure attempts to minimize the difference in short oligo frequency (summed over all oligos) between target and background data sets.  It calculates the desired weights for each background sequence to help minimize the error.  Due to the complexity of the problem, HOMER uses a simple hill-climbing approach by making small adjustment in background weight at a time.  It also penalizes large changes in background weight to avoid trivial solutions that a increase or decrease the weights of outlier sequences to extreme values.  The length of short oligos is controlled by the "-nlen <#>" option.

