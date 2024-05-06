
# full_stack_ChromHMM_annotations
Data of genome annotation from full-stack ChromHMM model trained with 1032 datasets from 127 reference epigenomes. Please refer to the published manuscript at <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02572-z"> Genome Biology </a>
# Download links:
Data of full-stack genome annotations for reference assemblies hg19 can be found <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//full_stack/full_stack_annotation_public_release/hg38"> here </a>. Within this folder:
```
├── hg19_genome_100_browser.bb: big bed file of the hg19's full stack annotation. 
This file can only be viewed properly on the ucsc genome browser through the trackhub link:
https://public.hoffman2.idre.ucla.edu/ernst/2K9RS///full_stack/full_stack_annotation_public_release/hub.txt
├── hg19_genome_100_browser.bed.gz: bed.gz version of the bigbed file hg19_genome_100_browser.bb. 
This file can be viewed as a custom track on UCSC genome browser
├── hg19_genome_100_segments.bed.gz: bed.gz version of the segmentation with only 4 columns corresponding to 
chrom, start, end, full_stack state
├── state_annotation_processed_publish.csv: annotations of the states 
(a more complete, excel-format version of this file is Additional File 3 in our published paper) 
├── state_annot_README.md: read me file for file state_annotation_processed_publish.csv
└── trackDb.txt
```
Data of full-stack genome annotations for reference assemblies hg38 can be found <a href="https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//full_stack/full_stack_annotation_public_release/hg38"> here </a>. Full-stack annotations in hg38 were created by lifting-over the annotation from hg10 to hg38. The whole liftOver pipeline's code is provide </a href="https://github.com/ernstlab/full_stack_ChromHMM_annotations/tree/main/chromHMM_utilities/liftOver"> here </a>. Within this folder:
```
├── hg38_genome_100_browser.bb: bigbed file of the hg38's full-stack annotation lifted over from hg19. 
This file can only be viewed properly on the ucsc genome browser through the trackhub link:
https://public.hoffman2.idre.ucla.edu/ernst/2K9RS///full_stack/full_stack_annotation_public_release/hub.txt
├── hg38_genome_100_browser.bed.gz: bed.gz version fo the bigbed file hg38_genome_100_browser.bb. 
This file can be viewed as a custom track on UCSC genome browser
├── hg38_genome_100_segments.bed.gz: bed.gz version of the segmentation with only 4 columns corresponding to 
chrom, start, end, full_stack state
├── state_annotations_processed.csv: annotations of the states 
(a more complete, excel-format version of this file is Additional File 3 in our published paper)
├── state_annot_README.md: read me file for file state_annotation_processed_publish.csv
└── trackDb.txt
```
- Detailed description of states can be found at Additional File 3 in the  </a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02572-z"> manuscript </a>. If you want a csv file outline similar information about the states' descriptions, you can also use file </a href="https://public.hoffman2.idre.ucla.edu/ernst/UUKP7/state_annotations_processed.csv">state_annotations_processed.csv</a> provided in this folder.
- An analogous annotation in mouse can be found on </a href="https://github.com/ernstlab/mouse_fullStack_annotations">this page</a>.
# Notes about the versions of full-stack annotation:
Since the publication of our data, the annotation of full-stack annotation in **hg19** have **NOT** been changed (__except for our state names being changed from 0-based indexing system 0_GapArtf1-99_TSS2 to 1_GapArtf1-100_TSS2__). The **hg38** annotation have been changed in 3 version**: 
- Version 1 (from publication of this repository to August 11 2023): In this version, we tried to liftOver the annotation from hg19 to hg38 such that no bp in hg38 were mapped from multiple bps in hg19. However, it was brought to our attention in February 2024 that there are still some regions in hg38 that were mapped from multiple regions in hg19 in this version of the annotation. This is caused by an implicit assumption in our code that each 200-bp segment from hg19 will be mapped to another segment of length 200-bp in hg38 by the liftOver tool provided by UCSC genome browser. This is not the case, so in version 3 (May 6 2024 onwards), we fixed this issue with updated code in the </a href="https://github.com/ernstlab/full_stack_ChromHMM_annotations/tree/main/chromHMM_utilities/liftOver"> liftOver pipeline</a>. In this first version, the states were numbered by the 0-based system (0_GapArtf1-99_TSS2, but note that the names GapArtf1, TSS2, etc. are exactly how they appeared in the paper). In the current version, we fixed this issue and eliminated 2,979,514 bp in verion1’s annotation that were wrongly mapped from different regions in hg19 (0.09% of the hg38 assembly genome). 
- Version 3 (May 6 2024 onwards): In this version, the state annotation for hg38 is corrected such that there are absolutely no bps in hg38 that were mapped from multiple bps in hg19. The states are numbered by the 1-based system (1_GapArtf1-100_TSS2, but note that the names GapArtf1, TSS2, etc. are exactly how they appeared in the paper). The change from 0-based to 1-based state numbering was made for consistency with other universal chromatin state model annotations. 
- Version 2: (From August 11 2023___ to February 26 2024___ ) : We mistakenly shared the links for the data in this version such that the annotations hg38 are raw output of the UCSC Genome browser's liftOver software (missing the step of filtering out bps in hg38 that were mapped from multiple bps in hg19). 
As of May 6 2024, we have archived the data of version3 onto https://zenodo.org/records/5759119
# UCSC Genome browser track hub link:
You can view the full-stack annotations in hg19, hg38 and mm10 (see </a href="https://github.com/ernstlab/mouse_fullStack_annotations"> Mouse's Github </a>) as a trackhub on UCSC genome browser using the link: https://public.hoffman2.idre.ucla.edu/ernst/2K9RS///full_stack/full_stack_annotation_public_release/hub.txt
# Folders:
Within each subfolders inside this folder, there are readme that can help you understand and apply the code. Note: AF stands for Additional File
- chromHMM_utilities contains subfolders and files that are helpful in processing the output of chromHMM (processing the emission matrix, overlap enrichment, neighborhood enrichment, calculating the AUROCs to measure how an annotation can help recover the genomic location of a genome element of interest). 
- enrichment_with_roadmap_151825_states_analysis: contains code to reproduce the results presented in AF1: Fig. S8 and excel spreadsheet AF6
- gene_ontology_analysis contains code to reproduce the results of gene ontology enrichment analysis for all the full-stack states, as presented in AF2.
- replicate_ernst_etal_2011_supp_fig2: contains code to reproduce the results presented in AF1: Fig. S10-11 showing the average gene expression for each full-stack state.
- test_celltype_specificity: contains code to reproduce the results presented in AF1: Fig. S7 showing the resulting of statistical tests on the cell-type specificity of full-stack states, based on emission parameters.
- random_represent_fullStack_state_with_25states: contains code to reproduce results presented in AF1: Fig. S9 showing the estimated probabilities of overlap between a full stack state and a 25-state from concatenated annotations for each of the biosamples. 
# License:
All code is provided under the MIT Open Acess License
Copyright 2021 Ha Vu and Jason Ernst
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# Contact:
We tried to commented the code and provide as much details on how to reproduce the results as possible. If you run into problems, please contact Ha Vu (havu73@ucla.edu) 

