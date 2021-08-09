import pandas as pd
state_type_color_dict = {'quescient': 'white', 'HET': 'light puple', 'enhancers': 'orange', 'znf': 'aquamarine', 'promoters': '#ff4500', 'transcription' : '#006400', 'others': 'seashell', 'weak enhancers': 'yellow', 'polycomb repressed': 'silver', 'transcribed and enhancer': '#ADFF2F', 'DNase': '#fff44f', 'exon': '#3cb371', 'acetylations': 'lemonchiffon', 'weak promoters': 'purple', 'TSS': 'red', 'weak transcription' : '#228B22'}

state_annotation_fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/state_annotations.csv'

def get_state_color(row_data):
	return state_type_color_dict[row_data['state_type']]  # return the color that correspond to the state type
df = pd.read_csv(state_annotation_fn, sep = ',', header = 0)
df['color'] = df.apply(get_state_color, axis = 1) # apply function to each row
df.to_csv(state_annotation_fn, header = True, index = False, sep = ',')
print "Done"