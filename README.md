# GeneralNormalizer
# UNDER CONSTRUCTION
R set of functions to normalise any table.
Shortly these functions will be provided as an R Library associated with a Publication under revision.
To use the functions it is possible to create a conda environment like:

conda create -n "env name" r-gridextra r-wesanderson r-dplyr r-data.table r-fitdistrplus r-ggplot2 r-mixsmsn r-factominer bioconductor-edger bioconductor-deseq2 bioconductor-complexheatmap bioconductor-biocparallel

This will include all the key packages to run the sets of functions

The main Function is called: "RunNorm" which takes in a data table and a design table and perform in parallel normalisation via Skewed-Normal mixed distribution.

A vignette and case studies will be added soon.

For now the function has been tested for bulk RNAseq, ATACseq (and DNAse-seq), Chipseq of histone modifications (not TFs), Proteomics, degraded proteomics and Metabolomics.

The main priciple stands out from the observation that not always median normalisation is accurate as it assums that the probability to observe increases and decreases when comparing two samples is equal. We have observed that for a plethora of conditions and cellular states this assumpion does not hold true. 
This in combination with the observation that for particular conditions we might have that the majority of the observations will change when comparing different conditions.

Therefore we came up with a strategy to empirivcally identify invariant observations on the besis of the frequency distributions of their changes when comparing pair-wise samples. Those observation that will lie under a bellshaped cure (with different degree of symmetry and skewness after data normalisation) will be most likely invariant within the comparison as araising from the sum of different random effects (each following different distributions).
