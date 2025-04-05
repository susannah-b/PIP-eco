from pipeco import *
import pipeco
import re

class pathotype_alignment :
	def __init__( self ) :
		global path_data, path_software
		global dir_result, path_genome, dir_blast, dir_out
		path_data = os.path.join(os.getcwd(), "data") # Allows scripts to be run independently
		path_software = os.path.join(os.getcwd(), "software")
		path_output = os.path.join(os.getcwd(), "pipeco_out")
		path_genome0 = os.path.join(os.getcwd(), "input", "fasta")
		path_genome = path_genome0.replace( "/fasta", "/faa" )
		dir_out = path_output + "/01.marker_out"
		os.system( "mkdir -p %s" %dir_out )
	@staticmethod

	def marker_pathotype( ) :
		def marker_alignment_run( file_faa ) :
			os.system( "mkdir -p %s/alignment" %dir_out )
			file_udb = "%s/marker/marker_gene.udb" %path_data
			file_usearch = "%s/usearch" %path_software
			blast_out0 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".aln" ) )
			blast_out1 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".b6" ) )
			blast_out2 = "%s/alignment/%s" %( dir_out, file_faa.split( "/" )[ -1 ].replace( ".faa", ".m8" ) )
			m8_userfields = "-userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl"
			command = "%s -ublast %s -db %s -evalue 1e-9 -alnout %s -blast6out %s -userout %s %s" %( file_usearch, file_faa, file_udb, blast_out0, blast_out1, blast_out2, m8_userfields )
			#print ( command )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
			#print ( outtext )
		def align_parse1( grep_m8 ) :
			# path_temp = "/".join( grep_m8.split( "/" )[ :-1 ] )
			# w_cov_files = glob.glob( os.path.join( path_temp, "*_w_cov.m8" ) )
			# if not w_cov_files :
			grep_m8_file0 = grep_m8.split( "/" )[ -1 ].replace( ".m8", "_w_cov.m9" )
			grep_m8_file1 = grep_m8.split( "/" )[ -1 ].replace( ".m8", "_fil.tsv" )
			df_m8_edit0 = "%s/alignment/%s" %( dir_out, grep_m8_file0 )
			df_m8_edit1 = "%s/alignment/%s" %( dir_out, grep_m8_file1 )
			if os.path.exists( df_m8_edit0 ) and os.path.exists( df_m8_edit1 ) :
				print ( f"Files already exist: {df_m8_edit0} and {df_m8_edit1}" )
				return
			try :
				df_m8 = pd.read_csv( grep_m8, sep = "\t", header = None )
				#display( df_m8 )
				df_m8.columns = [ "query_locus", "target_hit", "pidentity", "length", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore", "qlen", "tlen" ]
				df_m8[ "coverage" ] = ( df_m8[ "length" ] / df_m8[ "tlen" ] ) * 100 
				df_m8 = df_m8.reset_index()
				df_m8[ "query_locus_drop" ] = df_m8.query_locus.str.split( " #" ).str[ 0 ]
				df_m8.loc[ df_m8.coverage > 100, "coverage" ] = 100.0
				df_m8_drop = df_m8[ [ "query_locus_drop", "target_hit", "bitscore", "evalue", "pidentity", "coverage" ] ]
				df_filter_group0 = df_m8_drop[ df_m8_drop.apply( lambda x: x.coverage >= 90, axis = 1 ) ]
				df_filter_group1 = df_m8_drop[ (df_m8_drop.coverage >= 50) & (df_m8_drop.pidentity >= 80) ]
				df_filter_group1 = df_filter_group1.reset_index( drop = True )
				df_fil_drop = df_filter_group1.drop_duplicates( [ "query_locus_drop" ], keep = "first" )
				df_filter_group0.to_csv( df_m8_edit0, index = False, sep = "\t" )
				df_fil_drop.to_csv( df_m8_edit1, index = False, sep = "\t", header = True )
			except Exception as e:
				print(f"Error in align_parse1: {e}")
				shutil.copy( grep_m8, df_m8_edit0 )
				shutil.copy( grep_m8, df_m8_edit1 )

		def align_parse2( grep_tsv ):
			grep_tsv_file1 = grep_tsv.split( "/" )[ -1 ].replace( "_fil.tsv", "_fil_pid.tsv" )
			df_tsv_edit1 = "%s/alignment/%s" %( dir_out, grep_tsv_file1 )
			#df_m8_edit0 = "%s/%s/blastp_results/%s" %( path_out0, comparision0, grep_m8_file0 )
			#df_m8_edit1 = "%s/%s/blastp_results/%s" %( path_out0, comparision0, grep_m8_file1 )
			if os.path.exists( df_tsv_edit1 ) :
				print ( f"Files already exist: {df_tsv_edit1}" )
				return
			try :
				df_tsv = pd.read_csv( grep_tsv, sep = "\t", header = None )
				df_tsv.columns = [ "query_locus_drop", "target_hit", "bitscore", "evalue", "pidentity", "coverage" ]
				#df_m8.columns = [ "query_locus", "target_hit", "pidentity", "length", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore", "qlen", "tlen" ]
				df_tsv['pidentity'] = pd.to_numeric(df_tsv['pidentity'], errors='coerce')
				df_filter_group = df_tsv[ df_tsv.apply( lambda x: x.pidentity >= 80, axis = 1 ) ]
				df_filter_group = df_filter_group.reset_index( drop = True )
				df_fil_drop = df_filter_group.drop_duplicates( [ "query_locus_drop" ], keep = "first" )
				#df_m8_drop = df_m8[ [ "query_locus_drop", "target_hit", "bitscore", "evalue", "pidentity", "coverage" ] ]
				df_fil_drop.to_csv( df_tsv_edit1, index = False, sep = "\t" )
			except : 
				shutil.copy( grep_tsv, df_tsv_edit1 ) 

		def categorize_pathotype( dic_hit_marker ) :
			def check_bfp( value ) :
				bfp_markers = ["10_marker", "11_marker", "12_marker", "13_marker", "14_marker", "15_marker", "16_marker", "17_marker", "18_marker", "19_marker", "20_marker", "21_marker", "22_marker", "23_marker"]
				return sum( value.get( marker, 0 ) != 0 for marker in bfp_markers ) >= 2
			def check_eae( value ) :
				eae_markers = ["24_marker", "25_marker", "26_marker", "27_marker", "28_marker", "29_marker"]
				return any( value.get( marker, 0 ) != 0 for marker in eae_markers )
			def check_stx( value ) :
				stx_markers = ["41_marker", "42_marker", "43_marker", "44_marker"]
				return any( value.get( marker, 0 ) != 0 for marker in stx_markers )
			def check_eaec( value ) :
				eaec_markers = ["1_marker", "2_marker", "3_marker", "4_marker"]
				return any( value.get( marker, 0 ) != 0 for marker in eaec_markers )
			def check_daec( value ) :
				daec_markers = ["5_marker", "6_marker", "7_marker", "8_marker", "9_marker"]
				return any( value.get( marker, 0 ) != 0 for marker in daec_markers )
			def check_eiec( value ) :
				eiec_markers = ["30_marker", "31_marker", "32_marker", "33_marker", "34_marker", "35_marker"]
				return any( value.get( marker, 0 ) != 0 for marker in eiec_markers )
			def check_etec( value ) :
				etec_markers = ["36_marker", "37_marker", "38_marker", "39_marker", "40_marker"]
				return any( value.get( marker, 0 ) != 0 for marker in etec_markers )
			def check_iro( value ) :
				iro_markers = ["48_marker", "49_marker", "50_marker", "51_marker"]
				return sum( value.get( marker, 0 ) != 0 for marker in iro_markers ) >= 2
			def check_expec_aiec( value ) :
				expec_aiec_markers = ["45_marker", "46_marker", "47_marker", "52_marker", "53_marker", "54_marker", "55_marker"]
				return check_iro( value ) or any( value.get( marker, 0 ) != 0 for marker in expec_aiec_markers )
			category_counts = {"EAEC": [], "EIEC": [], "DAEC": [], "EPEC": [], "ETEC": [], "EHEC": [], "STEC": [], "ExPEC_AIEC": [], "hybrid": {}, "none": []} # Removed duplicate
			for key, value in dic_hit_marker.items() :
				if all( val == 0 for val in value.values() ) :
					category_counts[ "none" ].append( key )
					continue
				matched_categories = list()
				if check_eae( value ) :
					if check_stx( value ) :
						matched_categories.append( "EHEC" )
					else :
						if check_bfp( value ) :
							matched_categories.append( "EPEC" )
						else :
							matched_categories.append( "EPEC" )
				elif check_stx( value ) :
					matched_categories.append( "STEC" )
				if check_eaec( value ) :
					matched_categories.append( "EAEC" )
				if check_daec( value ) :
					matched_categories.append( "DAEC" )
				if check_eiec( value ) :
					matched_categories.append( "EIEC" )
				if check_etec( value ) :
					matched_categories.append( "ETEC" )
				if check_expec_aiec( value ) :
					matched_categories.append( "ExPEC_AIEC" )
				if len( matched_categories ) == 1 :
					category_counts[ matched_categories[ 0 ] ].append( key )
				elif len( matched_categories ) > 1 :
					category_counts[ "hybrid" ][ key ] = matched_categories
				else :
					category_counts[ "none" ].append( key )
			return category_counts
			
		num_cores = multiprocessing.cpu_count()
		lst_genome = glob.glob( path_genome + "/*.faa" )
		for genome0 in tqdm( lst_genome, desc = "Processing MARKER GENE: alignment" ) :		
			marker_alignment_run( genome0 )
		lst_m8 = glob.glob( dir_out + "/alignment/*.m8" )
		for m8_file in tqdm( lst_m8, desc = "Processing MARKER GENE: pre-parse_1" ) :
			align_parse1( m8_file )
		lst_tsv = glob.glob( dir_out + "/alignment/*_fil.tsv" )
		for tsv_file in tqdm( lst_tsv, desc = "Processing MARKER GENE: pre-parse_2" ) :
			align_parse2( tsv_file )
		
		lst_tsv0 = glob.glob( dir_out + "/alignment/*fil_pid.tsv" )
		markers = { f"{i}_marker":0 for i in range( 1, 56 ) }
		dic_hit_marker = dict()
		for tsv in lst_tsv0 :
			strain = tsv.split( "/" )[ -1 ].replace( "_fil_pid.tsv", "" )
			dic_hit_marker[ strain ] = markers.copy()
			df_tsv = pd.read_csv( tsv, sep = "\t" )
			#display( df_tsv )
			for idx, row in df_tsv.iterrows() :
				hit = row[ "target_hit" ].split( "|" )[ 0 ]; pidentity = row[ "pidentity" ]
				if dic_hit_marker[ strain ][ hit ] == 0 :
					dic_hit_marker[ strain ][ hit ] = [ pidentity ]
				else : dic_hit_marker[ strain ][ hit ].append( pidentity )
			#break
		for strain, markers in dic_hit_marker.items() :
			for marker, pidentities in markers.items() :
				if isinstance( pidentities, list) and pidentities :
					avg_pidentity = np.mean( pidentities )
					dic_hit_marker[ strain ][ marker ] = avg_pidentity
		
		category_counts = categorize_pathotype( dic_hit_marker )
		df_marker = pd.DataFrame( dic_hit_marker ).T
		
		# Note: samples expected as 'GCF_[Number]' format otherwise this will fail
		df_marker[ "index_num" ] = df_marker.index.str.extract(r'(\d+)', expand = False ).astype( int )
		df_marker = df_marker.sort_values( by = "index_num", ascending = True ).drop(columns = "index_num" )
		df_marker.to_csv( dir_out + "/01.marker_out.tsv", sep = "\t" )
		category_counts_json = json.dumps( category_counts, indent = 4 )
		with open( dir_out + "/01.marker_out.json", "w") as file:
			file.write( category_counts_json )
		scatter_data = df_marker.reset_index().melt( id_vars = "index", var_name = "marker", value_name = "value" )
		scatter_data.rename( columns={"index": "strain"}, inplace = True )
		scatter_data[ "strain_num" ] = scatter_data[ "strain" ].str.extract(r'(\d+)', expand = False ).astype( int )
		scatter_data[ "marker_num" ] = scatter_data[ "marker" ].str.extract(r'(\d+)', expand = False ).astype( int )

		scatter_data = scatter_data.sort_values( by = [ "marker_num", "strain_num" ] )

		# Modified graph: group by pathotype and adjust spacing
		# Define the desired pathotype order
		pathotype_order = ["EAEC", "EIEC", "DAEC", "EPEC", "ETEC", "EHEC", "STEC", "ExPEC_AIEC", "hybrid", "none"]
		# Generate sorted_strains grouped by pathotype
		sorted_strains = []
		for pathotype in pathotype_order:
			if pathotype == "hybrid":
				strains = list(category_counts["hybrid"].keys())
			else:
				strains = category_counts[pathotype]
			# Sort strains within the pathotype by numeric value
			strains_sorted = sorted(strains, key=lambda x: int(re.search(r'G.*_(\d+)', x).group(1)))
			sorted_strains.extend(strains_sorted)
		# Determine start/end indices for each pathotype group
		group_positions = {}
		current_idx = 0
		for pathotype in pathotype_order:
			if pathotype == "hybrid":
				count = len(category_counts["hybrid"])
			else:
				count = len(category_counts[pathotype])
			if count == 0:
				continue
			start = current_idx
			end = current_idx + count - 1
			group_positions[pathotype] = (start, end)
			current_idx += count

		# Create a categorical ordering based on the sorted strain names
		strain_order = pd.CategoricalDtype(categories=sorted_strains, ordered=True)
		scatter_data["strain_cat"] = scatter_data["strain"].astype(strain_order)

		zero_data = scatter_data[ scatter_data[ "value" ] == 0 ]
		non_zero_data = scatter_data[ scatter_data[ "value" ] > 0 ]
		colors = [ "#3679F5", "#55D557", "#FF0404" ]
		n_bins = 100; cmap_name = 'custom_cmap'
		cm = LinearSegmentedColormap.from_list( cmap_name, colors, N = n_bins )

		plt.figure( figsize = ( 30, 25 ) )
		sns.scatterplot( data = zero_data, y = "marker_num", x = "strain_cat", marker = 'x', color = 'black', s = 20 )
		sns.scatterplot( data = non_zero_data, y = "marker_num", x = "strain_cat", size = "value", hue = "value", palette = cm, sizes = (100, 200) )
		plt.yticks( ticks = scatter_data[ "marker_num" ].unique(), labels = scatter_data[ "marker" ].unique() )
		plt.xticks(rotation=45, ha='right', fontsize=10)
		for i, strain in enumerate(sorted_strains):
			plt.axvline(i, color='grey', linestyle=':', linewidth=0.5, alpha=0.5)

		# Add brackets and labels
		bracket_y = -9
		label_y = -9.5
		for pathotype, (start, end) in group_positions.items():
			plt.annotate('', xy=(start, bracket_y), xytext=(end, bracket_y), arrowprops=dict(arrowstyle='|-|,widthA=0.5,widthB=0.5',
																				color='black', lw=1.5), annotation_clip=False)
			plt.annotate(pathotype, xy=((start + end)/2, label_y), ha='center', va='top', annotation_clip=False)
		plt.subplots_adjust(bottom=0.25)

		plt.title( "Marker gene alignment" )
		plt.ylabel( "Marker" , fontsize=18)
		plt.xlabel( "Strain" , fontsize=18)
		plt.legend( title = "Value" )
		plt.savefig( dir_out + "/marker_alignment_result.svg" )
		plt.savefig( dir_out + "/marker_alignment_result.png" )

if __name__ == "__main__":
    pa = pathotype_alignment()
    pa.marker_pathotype()