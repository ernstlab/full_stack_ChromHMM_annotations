# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# parameters about the number of marks
TOTAL_NUM_MARKS <- 1032
# metadata about marks' cell type
metadata_fn <- '../ROADMAP_metadata_july2013.csv'
metadata <- as.data.frame(read.csv(metadata_fn, sep = ","))
metadata <- metadata %>% select("Epigenome.ID..EID.", "GROUP", "COLOR", "TYPE", "Epigenome.name..from.EDACC.Release.9.directory.", "ANATOMY", "TYPE") %>% rename("CT_NAME" = "Epigenome.ID..EID.") %>% rename("Epig_name" = "Epigenome.name..from.EDACC.Release.9.directory.")
metadata <- metadata[!apply(is.na(metadata) | metadata == "", 1, all),]
# metadata has CT_NAME, GROUP, COLOR
all_cell_groups <- as.character(unique(metadata$metada.GROUP))
get_mark_ct <- function(mark_name){
  return (unlist(strsplit(mark_name, "-"))[1])
}

get_chrom_mark_name <- function(mark_name){
  chrom_mark <- unlist(strsplit(mark_name, "-"))[2]
  acetyl_mark_list <- c('H2AK9ac', 'H4K12ac', 'H2BK20ac', 'H3K56ac', 'H4K5ac', 'K3BK15ac', 'H2BK12ac', 'H3K14ac', 'H4K91ac', 'H2AK5ac', 'H2BK120ac', 'H2BK5ac', 'H3K18ac', 'H3K23ac', 'H3K4ac', 'H2BK15ac', 'H4K8ac')
  other_mark_list <- c('H3T11ph', 'H3K23me2', 'H3K9me1')
  if (chrom_mark %in% acetyl_mark_list){
        chrom_mark <- 'class1_acetyl'
  }
  if (chrom_mark %in% other_mark_list){
        chrom_mark <- 'others'
  }
  return(chrom_mark)
}

paste_col_name <- function(num, prefix){
  return(paste0(prefix, num))
}

# Function to get rid of the "_mean" in "S01_mean"
fix_colname_mean_emission_df <- function(one_col_name){
        state_name <- unlist(strsplit(one_col_name, "_"))[1]
        return(as.integer(substr(state_name, 2, nchar(state_name))))
}

get_state_name_right <- function(state_index){
        return(paste0("S", state_index))
}