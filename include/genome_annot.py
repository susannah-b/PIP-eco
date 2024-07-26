from pipeco import *
import pipeco

class gannot :
	def __init__( self ) :
		pass
	@staticmethod

	def prokka_sc( ) :
		def file_move( lst_file, path_from, path_to ) :
			for file in lst_file :
				src_file0 = "%s/%s" %( path_from, file )
				des_file0 = "%s/%s" %( path_to, file )
				shutil.copy( src_file0, des_file0 )
		num_cores = multiprocessing.cpu_count()
		#print ( num_cores )
		path_genome = pipeco.path_input
		lst_genome = os.listdir( path_genome )
		for genome0 in tqdm( lst_genome, desc = "Processing genome annotation" ) :
			genome = "%s/%s" %( path_genome, genome0 )
			if ".fna" in genome :
				genome_name = genome0.replace( ".fna", "" )
			if ".fasta" in genome :
				genome_name = genome0.replace( ".fasta", "" )
			genome_name = "_".join( genome0.replace( ".fna", "" ).split( "_" )[ 0:2 ] )
			path_temp_out0 = path_genome.replace( "fasta", "temp_anno" )
			path_out0 = path_genome.replace( "fasta", "gff" )
			path_out1 = path_genome.replace( "fasta", "faa" )
			path_out2 = path_genome.replace( "fasta", "ffn" )
			os.system( "mkdir -p %s" %path_temp_out0 )
			os.system( "mkdir -p %s" %path_out0 )
			os.system( "mkdir -p %s" %path_out1 )
			os.system( "mkdir -p %s" %path_out2 )
			prokka_command = "prokka --outdir %s --prefix %s --force --cpus 32 %s" %( path_temp_out0, genome_name, genome )
			#print ( prokka_command )
			result = subprocess.run( prokka_command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True )
			if result.returncode != 0 :
				print( f"Error occurred while running prokka for {genome0}" )
				print( result.stderr.decode() )
			lst_gff = [ gff_file for gff_file in os.listdir( path_temp_out0 ) if gff_file.endswith( ".gff" ) ]
			lst_faa = [ faa_file for faa_file in os.listdir( path_temp_out0 ) if faa_file.endswith( ".faa" ) ]
			lst_ffn = [ faa_file for faa_file in os.listdir( path_temp_out0 ) if faa_file.endswith( ".ffn" ) ]
			file_move( lst_gff, path_temp_out0, path_out0 )
			file_move( lst_faa, path_temp_out0, path_out1 )
			file_move( lst_ffn, path_temp_out0, path_out2 )

			#print ( result.stderr )
		# dir_wgs = pipeco.dir_wgs
		# dir_genome = pipeco.dir_genome
		# file_wgs = dir_wgs + "/*.fna"
		# lst_fna = glob.glob( file_wgs ) 
		# lst_fna.sort()
		# if lst_fna == [] :
		# 	print ( "\n**********\nFatal Error : No Genome sequence File ...\n on %s\n**********\n" % dir_wgs )
		# count_fna = 0
		# for file_wgs in lst_fna :
		# 	count_fna += 1
		# 	fna_id = "_".join( file_wgs.split( "/" )[ -1 ].replace(".fna", "" ).split( "_" )[ 0 : 2 ] )
		# 	print ( "Proceesing... %s" %fna_id )
		# 	path_out = "%s/%s" %( dir_genome, fna_id ) 
		# 	os.system( "mkdir -p %s" %path_out )

		# 	if count_fna % 10 == 0 :
		# 		clear_output()
		# 	command = "prokka --outdir %s --prefix %s --force --cpus %s %s" %( path_out, fna_id, num_cores, file_wgs )
		# 	print ( command )
		# 	(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
		# 	print ( outtext )
