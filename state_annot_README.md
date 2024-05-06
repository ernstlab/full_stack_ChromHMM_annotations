The file state_annotation_processed_publish.csv outlines the states and their characteristics. They are simply a subset of columns in Additional File 3, tab 'characterize_full_stack_state', published in our paper (paper link: https://link.springer.com/article/10.1186/s13059-021-02572-z, AF3 link: https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-021-02572-z/MediaObjects/13059_2021_2572_MOESM3_ESM.xlsx ).
We still encourage you to look at our Additional file 3-5 published in the paper to have a comprehensive set of results and characterizations about the 100 full-stack states. But file state_annotation_processed_publish.csv is still provided as a short version of state characterization for your convenience. 

Columns in the file: 
- state_presented_in_paper: state names as presented in the paper. Format: <state_group><number_within_group>, example: GapArtf1. This column should match the column "States  presented in paper" in AF3, tab "characterize_full_stack_state" in our paper. 

- state_presented_in_bedFile: state names as presented in the .bed.gz file for download. Format:  <state_number><state_group><number_within_group>. State numbers are 1-based, so that states are  1_GapArtf1 --> 100_TSS2. The state number added to the front of state name is for your convenience in locating the states' characteristics in all our results (all our figures, and state characterizations in AF3-5).

- long_annotations: OUr comments about what each state may imply biologically. They are based on our understanding of the results from analyses outlined in the paper, presented in all figures and supplementary figures. We tried to create detailed comments about the states' biological implications so that you have a good first impression about the states. However, we would still encourage you to look at the states of your interests in details by browsing through all the tabs in Additional Files 3-5 in our paper. 

Read the following long note if you have confusions after looking through this column: 
We note that sometimes, in this columns, you encourter comments about some states, such as: "24_ReprPC (H3K27me3) in Brain, Epithelial, muscles, Mesench,Neurosph, Mystat, IM90, Adipose. Others quescient" . This means that state ReprPC6 in our paper shows the most enrichments to the state 24_ReprPC (polycomb repressed element marked by only H3K27me3 presence) in a multiple cell types, and are most enriched with the quiescent state in other cell types. Here, the state 24_reprPC is actually one state in the per-cell-type model that were used in Roadmap to annotate 127 biosamples (https://egg2.wustl.edu/roadmap/web_portal/imputed.html#chr_imp). This comment is based on our enrichment analyses between the full-stack states with 25-state per-cell-type annotation, shown in Supp. Fig. 8-9 (AF1) and also in excel form in AF5. We analyzed, for each full-stack state, their enrichments and probabilities of overlapping each of the 25 per-cell-type states in the annotation of 127 biosamples from Roadmap. This analyses help us understand whether each full-stack state is associated with any chromatin states that are used to annotate individual cell type/ biosample. 

- group
- color
- notes_about_enr_with_repeats: our own notes about the states' enrichments with different repeat classes and families, as presented in Supplementary figure 28 and AF4.
- notes_about_enr_with_perCT_annot: our own notes about the states' enrichments with per-cell-type annotations, as presented in Supplementary figures 8-9 and AF5. 

