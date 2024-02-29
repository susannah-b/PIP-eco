from pipeco import *
import pipeco
import re

class pai_analysis0 :
	def __init__( self ) :
		global path_data, path_program
		global dir_result, dir_genome, dir_blast, dir_pai, dir_out
		path_data = pipeco.path_data
		path_program = pipeco.path_program
		dir_genome = pipeco.dir_genome
		dir_blast = pipeco.dir_blast
		dir_pai = pipeco.dir_pai
		dir_result = pipeco.dir_result
		dir_out = dir_result + "/pai_out"	
		os.system( "mkdir -p %s" %dir_out )
	@staticmethod

	def pai_align( ) :
		def usearch_run( file_usearch, file_db, out_dir, file_faa ) : 
			blast_out0 = "%s/%s" %( out_dir, file_faa.split( "/" )[ -1 ].replace( ".faa", ".aln" ) )
			blast_out1 = "%s/%s" %( out_dir, file_faa.split( "/" )[ -1 ].replace( ".faa", ".b6" ) )
			blast_out2 = "%s/%s" %( out_dir, file_faa.split( "/" )[ -1 ].replace( ".faa", ".m8" ) )

			m8_userfields = "-userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl"
			command = "%s -ublast %s -db %s -evalue 1e-9 -alnout %s -blast6out %s -userout %s %s" %( file_usearch, file_faa, file_db, blast_out0, blast_out1, blast_out2, m8_userfields )
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
				df_filter_group1 = df_m8_drop[ (df_m8_drop.coverage >= 90) & (df_m8_drop.pidentity >= 95) ]
				df_filter_group1 = df_filter_group1.reset_index( drop = True )
				df_fil_drop = df_filter_group1.drop_duplicates( [ "query_locus_drop" ], keep = "first" )
				df_filter_group0.to_csv( df_m8_edit0, index = False, sep = "\t" )
				df_fil_drop.to_csv( df_m8_edit1, index = False, sep = "\t" )
			except :
				shutil.copy( grep_m8, df_m8_edit0 )
				shutil.copy( grep_m8, df_m8_edit1 )
		def pai_landscape( in_file, out_file ) :
			dic_pai = { "AGI-1": 0, "PAI_II_APEC-O1": 0, "AGI-3": 0, "PAI-I_AL862": 0, "LEE": 0, "espCPAI": 0, "ETT2": 0, "GimA": 0, "HPI": 0, "OI-122": 0, "LEEII": 0, "LIM": 0, "LPA": 0, "Notname": 0, "PAI_III_536": 0, "PAI_I_APEC-O1": 0, "OI-57": 0, "PAI_I_4787": 0, "PAI_I_536": 0, "PAI_V_536": 0, "PAI_II_536": 0, "PAI_I_CFT073": 0, "PAI_I_CL3": 0, "PAI_II_4787": 0, "PAI_II_CFT073": 0, "PAI_III_APEC-O1": 0, "PAI_IV_536": 0, "PAI_IV_APEC-O1": 0, "TAI": 0, "PAI-I-AL862": 0, "SE-PAI": 0, "Tia-PAI": 0 }
			df_hit = pd.read_csv( in_file, sep = "\t" )
			for idx, row in df_hit.iterrows() :
				pai_hit = row[ "target_hit" ]; pid = row[ "pidentity" ]
				hit_pai = pai_hit.split( "|" )[ 3 ].split( "_[" )[ 0 ].replace( " ", "_" )
				dic_pai[ hit_pai ] += 1
			df = pd.DataFrame.from_dict( dic_pai, orient = "index", columns = [ "count" ] ).reset_index()
			df.rename( columns = { "index" : "pai" } , inplace = True ) 
			df_t = df.set_index( "pai" )
			log_df = df_t.applymap( lambda x: np.log10( x ) if x > 0 else 0 )
			fig, ( ax1, ax2 ) = plt.subplots( 1, 2, figsize = ( 10, 15 ), gridspec_kw={ "width_ratios": [ 1, 3 ] } )
			sns.heatmap( log_df, annot = True, linecolor = "black", linewidths = .1, fmt = ".2f", cmap = "Greys", cbar = False, ax = ax1 )
			barplot = sns.barplot( y = df[ "pai" ], x = df[ "count" ], ax = ax2, orient = "h" )
			for index, value in enumerate( df[ "count" ] ) :
				ax2.text( value, index, str( value ), color = "black", ha = "center", va = "center" )
			ax2.set_title( "PAI Counts" )
			ax2.set_xlabel( "Count" )
			plt.tight_layout()
			plt.savefig( out_file )
		
		num_cores = multiprocessing.cpu_count()
		file_usearch = "%s/usearch" %path_program
		file_db = "%s/UniqPAIdb/pai_determinant.udb" %path_data
		print ( num_cores )
		for genome in os.listdir( dir_genome ) :
			file_faa = "%s/%s/%s.faa" %( dir_genome, genome, genome )
			out_dir = "%s/%s" %( dir_pai, genome )
			out_dir1 = "%s/%s" %( dir_out, genome )
			os.system( "mkdir -p %s" %out_dir )
			os.system( "mkdir -p %s" %out_dir1 )
			usearch_run( file_usearch, file_db, out_dir, file_faa )
			file_m8 = "%s/%s.m8" %( out_dir, genome )
			usearch_parse( out_dir, file_m8 )
			file_tsv = "%s/%s_fil.tsv" %( out_dir, genome )
			out_fig = "%s/%s_pai_alignment.svg" %( out_dir1, genome )
			pai_landscape( file_tsv, out_fig )