from pipeco import *
import pipeco

class prokka :
	def __init__( self ) :
		pass
	@staticmethod

	def prokka_sc( ) :
		num_cores = multiprocessing.cpu_count()
		print ( num_cores )
		dir_wgs = pipeco.dir_wgs
		dir_genome = pipeco.dir_genome
		file_wgs = dir_wgs + "/*.fna"
		lst_fna = glob.glob( file_wgs ) 
		lst_fna.sort()
		if lst_fna == [] :
			print ( "\n**********\nFatal Error : No Genome sequence File ...\n on %s\n**********\n" % dir_wgs )
		count_fna = 0
		for file_wgs in lst_fna :
			count_fna += 1
			fna_id = "_".join( file_wgs.split( "/" )[ -1 ].replace(".fna", "" ).split( "_" )[ 0 : 2 ] )
			print ( "Proceesing... %s" %fna_id )
			path_out = "%s/%s" %( dir_genome, fna_id ) 
			os.system( "mkdir -p %s" %path_out )

			if count_fna % 10 == 0 :
				clear_output()
			command = "prokka --outdir %s --prefix %s --force --cpus %s %s" %( path_out, fna_id, num_cores, file_wgs )
			print ( command )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
			print ( outtext )
