# ECLIPSE

**ECLIPSE (Extraction of Clonality from LIquid bioPSiEs)** is an R package for clonal deconvolution of tumour-informed ctDNA data using clonality  and copy number information from tumour tissue. It allows calculation of changes in cancer cell fraction (growth and death) of subclones that have been found in tumour tissue samples. It also calculates measure of tumour burden in a patient which is independent of the copy number state of the targeted variants ( ie sample purity/cellularity ) and hence is comparable across different tumours which may have significant differences in DNA content per cell (for example if whole genome doubling has occurred). The package is designed and has been validated using high quality deep targeted multiplex-archored PCR sequencing of variants from multi-region exome sequencing of NCLSCs as part of the TRACERx project. In this dataset set ECLIPSE has been effective in samples with at least 0.1% purity, greatly extending the number of samples where by clonal deconvolution would be possible using deep exome-seq and informatic methods such as PyClone.

![](data-raw/ECLPSE_infographic.jpg)

## Installation & loading

You can use devtools::install_github() to install cloneMap from this repository:
Temporary authenication token for reviewers included

```R
devtools::install_github('amf71/ECLIPSE',auth_token = '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')`
```

load package:

```R
library(eclipse)
```

## Data inputs for ECLIPSE 

For each mutation ECLIPSE requires:
* mutation identifiers (Chromosome, position, reference allele, alternative allele)
* sample identifier
* number of supporting reads
* sequencing depth
* estimated background error rate
* clone identifier
* binary call for whether the mutation is clonal or subclonal
* mutation multiplicity
* total copy number at the mutated locus in tumour cells
* total copy number at the mutated locus in non-tumour cells (default = 2)  

ECLIPSE also takes several optional inputs:
* binary call for if a variant should be filtered for clone and tumour present calls due to high background error (the probability, for any given read, to observe the specified mutation due to sequencing error)
* binary call for if a variant  should be filtered from all analysis for a specific sample  
* measurement for the maximally expected normalised standard deviation of CCF in high confidence clones used to identify clones with incoherent CCF
distributions which may represent mutation clusters that are not true clones. This can be calculated using the `extract_normalised_sd` function which leverages higher tumour purity cfDNA samples in your data and clonal clusters to estimate this. 

The clone identifier, clonal vs subclonal status, mutation multiplicity and total copy number in tumour cells can be calculated using standard copy number extraction and clonal deconvolution methods (ASCAT, Battenberg, Pyclone, DpCLust) used for high tumour purity (>10%) samples, for example from tissue samples, which can then be used as estimates for these variables at the time of ctDNA sampling. Clonal status can be more accurately and comprehensively extracted from the sequencing of multiple high purity samples from the same patient, as is performed in TRACERx, but is not essential. See Application of ECLIPSE to the TRACERx cfDNA data section for further details. 


## Stepwise description of ECLIPSE 

1. **VAF denoising**. Variant allele frequencies (VAFs) are denoised by subtracting the estimated background error, provided to ECLIPSE for each variant. For a description of estimating background error in this dataset see MRD calling algorithm section. Variants in each clone are grouped into clusters (k-means clustering) with similar background error profiles where the number of groups is determined by the sum of the error estimated across all variants, so that if equally dividing the total error from all variants of a clone each group would have a combined error of at least one mutant read. Therefore, if a clone has a total combined error of less than two mutant reads only one cluster will be used. A maximum number of clusters is set to four as default which was used for this dataset. The average background error of each group per variant is subtracted from the number of supporting reads for all variants in each group and divided by the sequencing depth to calculate denoised VAFs. 
2. **ctDNA tumour purity calculation**. Deniosed VAFs are used with mutation multiplicities, total copy number at the mutated locus and clonal vs subclonal mutation status for each mutation provided to ECLIPSE to calculate an estimate of ctDNA tumour purity (see above) using the equation shown in Supplementary Figure 3b for each clonal mutation. The equation shown in sup figure 3c is a rearrangement of that shown in sup figure 3c for clonal mutations where CCF = 1. We summarise the mean of these values to provide a final estimate of ctDNA tumour purity per sample. 
3. **CCF calculation per mutation and subclone**.  For all mutations the sample’s ctDNA tumour purity, denoised VAF, multiplicity and total copy number at the mutated locus used in equation show in Supplementary Figure 3c to calculate an estimate of CCF for each mutation in a given plasma sample. The clone identities for each mutation are provided to ECLIPSE and should be calculated independently using standard methods which leverage SNP coverage applicable to high purity samples39–41. A mean of the per-mutation CCFs is used as a CCF estimate for each clone. 
4. **Poor quality clone identification**. Mutation clustering using standard methodologies are imperfect and are fitted to high purity samples (usually matched tissue samples) used for clustering, excluding lower purity samples which ECLIPSE is able to analyse with deep targeting sequencing. Clusters which are not real, i.e. for which the mutations are not in the same cells, may not continue to track at similar CCFs in data from new samples. To capture such clusters the distribution of ECLIPSE calculated CCFs in each clone in a ctDNA sample are quantified using normalised standard deviations (SDs). This can then be compared to the expected CCF distributions of high confidence clones, for example clonal clusters in higher purity plasma samples. In our data we quantified the normalised SD of all clonal clusters in samples of greater than 5% purity and took the upper 95% confidence interval for this data calculated at 0.56. Subclonal clusters with normalised SDs for CCFs > 0.56 were considered of poor quality and were not considered for analysis. This identified 2.6% of clones in the TRACERx data as of poor quality. Expected CCF distributions will be highly dependent on the input data for ECLIPSE and should therefore be benchmarked on each data set. A function in the ECLIPSE R package is provided to calculate an upper 95% CI of normalised SDs for CCFs in clonal clusters in high purity samples as was performed for this dataset. 
5. **Clone present calling**. Determination of whether each clone was present or absent from each sample (see High specificity subclone detection section) by comparing the sum of expected background error to the sum of observed signal across all variants in the subclone with a one-sided Poisson test. Mutations with high error that should be excluded from these calculations can be specified. 
6. **Tumour present calling**.  Determination of whether any tumour cells were present for each sample (see above) by comparing the summed expected background error to the summed observed signal across all variants tracked in the sample with a Poisson test. Mutations with high noise that should be excluded from these calculations can be specified. 
7. **Minimal detectable CCF estimation for each subclone**. Determination of the CCF equivalent to the minimal number of supporting reads across all variants in a subclone that would be required for a significant clone present Poisson test (P<0.01, see High specificity subclone detection section). 
8. **Minimal detectable CCF estimation for an average subclone for each sample**.  Determination of the CCF equivalent to the minimal number of supporting reads across all variants in a representative subclone that would be required for a significant clone present Poisson test (P<0.01, see High specificity subclone detection section). The background is taken as an average of the background error in all subclonal mutations tracked in a given sample and is representative for a subclone tracked by four mutations as default, the average number tracked in this dataset. This value allows comparisons of minimally detectable CCF limit across samples.
9. **Minimal detectable purity estimation for each sample**. Determination of the purity equivalent to the minimal number of supporting reads across all tracked variants that would be required for a significant tumour present Poisson test (P<0.01). 
10. **Testing for the absence of a complete clonal sweep for each subclone**. A subclone which is detected in high purity samples used for mutation clustering may expand through a full clonal sweep later in the disease course. We would therefore expect to observe CCFs of 100%, indistinguishable from CCFs of clonal mutations after such an event. For each subclone in each sample, a Wilcoxon test is performed to compare the CCFs of each subclone to the CCFs of clonal mutations in the same sample. The resulting P value indicates whether there is significant evidence that the subclone has significantly <100% CCF and therefore is only present in the minority of tumour cells, without a full clonal sweep. 








