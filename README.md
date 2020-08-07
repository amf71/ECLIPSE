# ECLIPSE

**ECLIPSE (Extraction of Clonality from LIquid bioPSiEs)** is an R package for clonal deconvolution of tumour-informed ctDNA data using clonality 
and copy number information from tumour tissue. It allows calculation of changes in cancer cell fraction (growth and death) of subclones that 
have been found in tumour tissue samples. It also calculates measure of tumour burden in a patient which is independent of the copy number 
state of the targeted variants ( tumour cellularity ) and hence is comparable across different tumours which may have significant differences 
in DNA content per cell (for example if whole genome doubling has occurred). The package is designed and has been validated using high quality 
deep amplicon sequencing of variants targeted from multi-region exome sequencing as part of the TRACERx project. In this dataset set it has 
been effective in samples with at least 0.1% tumour fraction greatly extending the number of samples where by clonal deconvolution would be 
possible using deep exome-seq and methods such as PyClone.



