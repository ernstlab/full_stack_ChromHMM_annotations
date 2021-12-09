
# full_stack_ChromHMM_annotations
Data of genome annotation from full-stack ChromHMM model trained with 1032 datasets from 127 reference epigenomes 
# Download links:
Data of full-stack genome annotations for reference assemblies hg19 and hg38 can be found <a href="https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/"> here</a>: 
Within this folder:
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/hg19_genome_100_segments.bed.gz">hg19_genome_100_segments.bed.gz</a> contains a simple four column .bed file of full-stack state annotation in hg19 assembly. Since our training data (1032 input data tracks) are in hg19, the assembly used for original training and annotation. The fourth column contains a state label with a prefix number that can be used to order the states. The OverlapEnrichment and NeighborhoodEnrichment commands of ChromHMM with the '-labels' option can compute enrichments for this file and order states based on the prefix number.
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/hg19_genome_100_browser.bed.gz">hg19_genome_100_browser.bed.gz</a> contains a browser file of full-stack state annotation in hg19 assembly. This file is compatible to for UCSC genome browser. Since our training data (1032 input data tracks) are in hg19, the assembly used for original training and annotation.
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/hg38lift_genome_100_segments.bed.gz">hg38lift_genome_100_segments.bed.gz</a> is similar to hg19_genome_100_segments.bed.gz, but based on a liftOver of each 200bp bin to hg38. Positions mapped to from multiple locations in hg19 did not receive an annotation. Since the 200-bp intervals are no longer maintained when using OverlapEnrichment and NeighborhoodEnrichment the options '-b' and '-lowmmem' should be used with this command in addition to the '-labels' option as with the hg19 file.
- File <a href="https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/hg38lift_genome_100_browser.bed.gz">hg38lift_genome_100_browser.bed.gz</a> is similar to hg19_genome_100_browser.bed.gz but is for the liftOver to hg38. Positions mapped to from multiple locations in hg19 did not receive an annotation.

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
