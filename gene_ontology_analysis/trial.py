import pandas as pd 

fn = '/u/home/h/havu73/project-ernst//data/housekeeping_genes/HK_genes.csv'
gdf = pd.read_csv(fn, header = None, index_col = None, sep = '\t')
gdf.columns = ['gene_name']
efn = '/u/home/h/havu73/project-ernst//data/roadmap_epigenome/gene_exp/trial'
edf = pd.read_csv(efn, header = None, index_col = None, sep = '\t')
edf = edf.dropna()
edf.columns = ['chrom', 'start', 'end', 'strand', 'gene_type', 'gene_name']
gdf['gene_name'] = gdf['gene_name'].apply(lambda x: str(x))
edf['gene_name'] = edf['gene_name'].apply(lambda x: str(x))
t = gdf.merge(edf, left_on = 'gene_name', right_on = 'gene_name')
t['chrom'] = t['chrom'].apply(lambda x: 'chr' + x)