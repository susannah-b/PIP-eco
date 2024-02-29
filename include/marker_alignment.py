from pipeco import *
import pipeco
import re

class local_alignment :
	def __init__( self ) :
		global path_data, path_program
		global dir_result, dir_genome, dir_blast, dir_out
		path_data = pipeco.path_data
		path_program = pipeco.path_program
		dir_genome = pipeco.dir_genome
		dir_blast = pipeco.dir_blast
		dir_result = pipeco.dir_result
		dir_out = dir_result + "/marker_out"	
		os.system( "mkdir -p %s" %dir_out )
	@staticmethod

	def blastp_run( ) :
		def usearch_run( file_usearch, file_marker_db, out_dir, file_faa ) :
			blast_out0 = "%s/%s" %( out_dir, file_faa.split( "/" )[ -1 ].replace( ".faa", ".aln" ) )
			blast_out1 = "%s/%s" %( out_dir, file_faa.split( "/" )[ -1 ].replace( ".faa", ".b6" ) )
			blast_out2 = "%s/%s" %( out_dir, file_faa.split( "/" )[ -1 ].replace( ".faa", ".m8" ) )

			m8_userfields = "-userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl"
			command = "%s -ublast %s -db %s -evalue 1e-9 -alnout %s -blast6out %s -userout %s %s" %( file_usearch, file_faa, file_marker_db, blast_out0, blast_out1, blast_out2, m8_userfields )
			print ( command )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
			print ( outtext )
		
		def usearch_parse( out_dir, grep_m8 ) :
			grep_m8_file0 = grep_m8.split( "/" )[ -1 ].replace( ".m8", "_w_cov.m8" )
			grep_m8_file1 = grep_m8.split( "/" )[ -1 ].replace( ".m8", "_fil.tsv" )
			df_m8_edit0 = "%s/%s" %( out_dir, grep_m8_file0 )
			df_m8_edit1 = "%s/%s" %( out_dir, grep_m8_file1 )
			try : 
				df_m8 = pd.read_csv( grep_m8, sep = "\t", header = None )
				df_m8.columns = [ "query_locus", "target_hit", "pidentity", "length", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore", "qlen", "tlen" ]
				df_m8[ "coverage" ] = ( df_m8[ "length" ] / df_m8[ "tlen" ] ) * 100 
				df_m8 = df_m8.reset_index()
				df_m8[ "query_locus_drop" ] = df_m8.query_locus.str.split( " #" ).str[ 0 ]
				df_m8.loc[ df_m8.coverage > 100, "coverage" ] = 100.0
				df_m8_drop = df_m8[ [ "query_locus_drop", "target_hit", "bitscore", "evalue", "pidentity", "coverage" ] ]
				df_filter_group0 = df_m8_drop[ df_m8_drop.apply( lambda x: x.coverage >= 90, axis = 1 ) ]
				df_filter_group1 = df_m8_drop[ (df_m8_drop.coverage >= 90) & (df_m8_drop.pidentity >= 80) ]
				df_filter_group1 = df_filter_group1.reset_index( drop = True )
				df_fil_drop = df_filter_group1.drop_duplicates( [ "query_locus_drop" ], keep = "first" )
				df_filter_group0.to_csv( df_m8_edit0, index = False, sep = "\t" )
				df_fil_drop.to_csv( df_m8_edit1, index = False, sep = "\t" )
			except :
				shutil.copy( grep_m8, df_m8_edit0 )
				shutil.copy( grep_m8, df_m8_edit1 )

		num_cores = multiprocessing.cpu_count()
		file_usearch = "%s/usearch" % path_program
		file_marker_db = "%s/marker_db/marker_gene_seq.udb" %path_data
		lst_genome = os.listdir( dir_genome ) 
		for genome0 in lst_genome :
			file_faa = "%s/%s/%s.faa" %( dir_genome, genome0, genome0 )
			out_dir = "%s/%s" %( dir_blast, genome0 )
			os.system( "mkdir -p %s" % out_dir )
			usearch_run( file_usearch, file_marker_db, out_dir, file_faa )
			file_m8 = "%s/%s.m8" %( out_dir, genome0 )
			usearch_parse( out_dir, file_m8 )

	@staticmethod
	def pathotype_assign( ) :
		def check_general_category( value, markers ) :
			return any( value.get( marker, 0 ) != 0 for marker in markers ) and all( value.get( marker, 0 ) == 0 for marker in value if marker not in markers )

		def check_ehec_stec( value ) :
			ehec_stec_markers_sets = [ 
			[ "11_marker", "15_marker", "27_marker" ], 
			[ "11_marker", "15_marker", "28_marker" ], 
			[ "11_marker", "27_marker" ],
			[ "11_marker", "28_marker" ], 
			[ "27_marker" ],
			[ "28_marker" ] ]
			return any( all( value.get( marker, 0 ) != 0 for marker in markers_set ) for markers_set in ehec_stec_markers_sets ) and all( value.get( marker, 0 ) == 0 for marker in value if marker not in [ "11_marker", "15_marker", "27_marker", "28_marker" ] )

		def classify_pathotype( row ) :
			pathotypes = list()
			if any( row[ f"{i}_marker" ] > 0 for i in range( 1, 6 ) ) :
				pathotypes.append( "EIEC" )
			if any( row[ f"{i}_marker" ] > 0 for i in range( 6, 11) ) :
				pathotypes.append( "ETEC" )
			if any( row[ f"{i}_marker" ] > 0 for i in [ 11, 12 ] ) :
				pathotypes.append( "EPEC" )
			if any( row[ f"{i}_marker" ] > 0 for i in list( range( 13, 15 ) ) + list( range( 16, 27 ) ) ) :
				pathotypes.append( "ExPEC_APEC" )
			if ( row["27_marker"] > 0 or row["28_marker"] > 0 ) :
				if "EPEC" in pathotypes :
					pathotypes.remove( "EPEC" )
				pathotypes.append( "STEC_EHEC" )
			if not pathotypes :
				return "Unknown" 
			return ", ".join( pathotypes )
		
		markers = { f"{i}_marker":0 for i in range( 1, 29 ) }
		dic_hit_marker = dict()
		for idx, genome in enumerate( os.listdir( dir_blast ) ) :
			file_align = glob.glob( "%s/%s/*.tsv" %( dir_blast, genome ) )
			dic_hit_marker[ genome ] = markers.copy()
			df_align = pd.read_csv( file_align[ 0 ], sep = "\t" )
			for i0, r0 in df_align.iterrows() :
				hit_marker = r0[ "target_hit" ].split( "|" )[ 0 ]
				hit_patho = r0[ "target_hit" ].split( "|" )[ 1 ]
				pidentity = r0[ "pidentity" ]
				if dic_hit_marker[ genome ][ hit_marker ] == 0 :
					dic_hit_marker[ genome ][ hit_marker ] = [ pidentity ]
				else : dic_hit_marker[ genome ][ hit_marker ].append( pidentity )
			idx += 1
		for strain, markers in dic_hit_marker.items() :
			for marker, pidentities in markers.items() :
				if isinstance( pidentities, list ) and pidentities :
					avg_pidentity = np.mean( pidentities )
					dic_hit_marker[ strain ][ marker ] = avg_pidentity
		#print( dic_hit_marker )
		marker_sets = {
		"EIEC": ["1_marker", "2_marker", "3_marker", "4_marker", "5_marker"],
		"ETEC": ["6_marker", "7_marker", "8_marker", "9_marker", "10_marker", "15_marker"],
		"EPEC": ["11_marker", "12_marker", "15_marker"],
		"ExPEC_APEC": ["13_marker", "14_marker", "15_marker", "16_marker", "17_marker", "18_marker", "19_marker", "20_marker", "21_marker", "22_marker", "23_marker", "24_marker", "25_marker", "26_marker"],
		"EHEC_STEC": ["11_marker", "15_marker", "27_marker", "28_marker"]
		}

		category_counts = { "EIEC": 0, "ETEC": 0, "EPEC": 0, "ExPEC_APEC": 0, "EHEC_STEC": 0, "hybrid": 0, "none": 0, "unknown" : 0 }
		category_keys = { category: [] for category in [ "EIEC", "ETEC", "EPEC", "ExPEC_APEC", "EHEC_STEC", "hybrid", "none", "unknown" ] }
		for key, value in dic_hit_marker.items() :
			if all( val == 0 for val in value.values() ) :
				category_counts[ "none" ] += 1
				category_keys[ "none" ].append( key )
				continue

			if all( value.get( marker, 0 ) == 0 for marker in value if marker != "15_marker" ) and value.get( "15_marker", 0 ) != 0 :
				category_counts[ "unknown" ] += 1
				category_keys[ "unknown" ].append( key )
				continue

			matched_categories = list()
			for category, markers in marker_sets.items() :
				if category == "EHEC_STEC" and check_ehec_stec( value ) :
					matched_categories.append( category )
				elif category != "EHEC_STEC" and check_general_category( value, markers ) :
					matched_categories.append( category )

			if len( matched_categories ) == 1 :
				category_counts[ matched_categories[ 0 ] ] += 1
				category_keys[ matched_categories[ 0 ] ].append( key )
			elif len( matched_categories ) > 1 :
				category_counts[ "hybrid" ] += 1
				category_keys[ "hybrid" ].append( key )
			else :
				category_counts[ "hybrid" ] += 1
				category_keys[ "hybrid" ].append( key )
		for category, keys in category_keys.items() :
			print ( f"{category}: {len(keys)} items" )
			print ( f"Genome: {keys}\n" )
		
		df_marker_out = pd.DataFrame( dic_hit_marker ).T
		#display( df_marker_out )
		df_marker_out.to_csv( dir_out + "/01.marker_out.tsv" , sep = "\t" )
		scatter_data = df_marker_out.reset_index().melt( id_vars = "index", var_name = "marker", value_name = "value" )
		scatter_data.rename( columns = { "index": "strain"}, inplace = True )
		scatter_data[ "marker_num" ] = scatter_data[ "marker" ].apply( lambda x: int( re.search(r'\d+', x).group() ) )
		scatter_data.sort_values( by = "marker_num", inplace = True )
		
		unique_strains = len( df_marker_out.index )
		height = 0.5
		total_height = unique_strains * height
		plt.figure(figsize=(15, max(5, total_height)))
		palette = sns.color_palette( "coolwarm", as_cmap = True )
		sns.scatterplot( data = scatter_data[ scatter_data[ "value" ] > 0 ],
		x = "marker", y = "strain", size = "value", sizes = ( 100, 500 ), 
		hue = "value", palette = palette, alpha = 0.9 )
		plt.scatter( scatter_data[ scatter_data[ "value" ] == 0 ][ "marker"],
			scatter_data[ scatter_data[ "value"] == 0 ][ "strain" ],
			marker = "x", s = 20, color = "gray" )
		x_labels = scatter_data.sort_values( by = "marker_num" )[ "marker" ].unique()
		plt.xticks( rotation = 45 )
		plt.ylim( -1, unique_strains )
		plt.tight_layout()
		plt.savefig( dir_out + "/01.marker_out.svg" )

		file_marker_out = dir_out + "/01.marker_out.tsv" 
		df_markers = pd.read_csv( file_marker_out, sep = "\t" )
		df_markers[ "Pathotype" ] = df_markers.apply( classify_pathotype, axis = 1 )
		pathotypes = [ "EIEC", "ETEC", "EPEC", "ExPEC_APEC", "STEC_EHEC" ]
		matrix = pd.DataFrame( 0, index = pathotypes, columns = pathotypes ) 
		display( df_markers )
		df_markers.to_csv( dir_out + "/01.marker_out.tsv", sep = "\t" )