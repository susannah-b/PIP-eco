from pipeco import *
import pipeco
import re

class pan_phylo :
	def __init__( self ) :
		global path_data, path_program
		global dir_result, dir_genome, dir_blast, dir_pan, dir_out
		path_data = pipeco.path_data
		path_program = pipeco.path_program
		dir_genome = pipeco.dir_genome
		dir_blast = pipeco.dir_blast
		dir_pan = pipeco.dir_pan
		dir_result = pipeco.dir_result
		dir_out = dir_result + "/phylo_out"	
		os.system( "mkdir -p %s" %dir_out )
	@staticmethod

	def phylo_analysis( ) :
		def transform_matrix( df ) :
			df_transposed = df.drop( "Gene", axis = 1 ).transpose()
			df_transposed = df_transposed.replace( 1, "A" )
			df_transposed = df_transposed.replace( 0, "T" )
			return df_transposed
		def create_matrix( df, file_out ) :
			with open( out_file, "w" ) as f :
				for strain, sequence in df.iterrows() :
					out_seq = "".join( sequence.astype( str ) )
					f.write( f">{strain}\n{out_seq}\n" )
			print( "Wrote... %s" %file_out )
		def panaroo_default( gff, dir_temp, num_cores ) :
			command = "panaroo -i %s -c 0.8 --len_dif_percent 1.0 -o %s -t %s --clean-mode moderate" %( gff, dir_temp, num_cores )
			print( command ) 
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command )
			print ( outtext )
		def panaroo_integrated( backbone, input_gff, dir_out ) :
			command0 = "panaroo-integrate -d %s -i %s -c 0.8 --len_dif_percent 1.0 -t 32 -o %s" %( backbone, input_genome, out_dir0 )
		def panaroo_merged( draft_merged_graph, backbone, dir_out, num_cores ) :
			command1 = "panaroo-merge -d %s %s -t %s -o %s" %( draft_merged_graph, backbone, num_cores, dir_out )
			print ( command1 )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command1 )
			print ( outtext )
		def msa_mafft( in_file, out_aln, num_cores ) :
			command2 = "mafft --auto --thread %s %s > %s" %( num_cores, in_file, out_aln )
			print ( command2 )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command2 )
			print( outtext )
		def make_phylo( file_muscle, out_aln, out_nwk ) :
			command3 = "%s -maketree -in %s -out %s -cluster neighborjoining" %( file_muscle, out_aln, out_nwk )
			print ( command3 )
			(exitstatus, outtext) = subprocess.getstatusoutput( "%s" % command3 )
			print ( outtext )
		def generate_tree( nwk, out_file ) : 
			with open( nwk, "r" ) as file : 
				data = file.read().replace( "\n", "" )
			canvas = toyplot.Canvas( width = 450, height = 850, style = { "background-color": "white" } )
			axes = canvas.cartesian()
			tre = toytree.tree( data )
			tre.draw( 
				axes = axes,
				tip_labels = True, 
				node_labels = True, 
				edge_type = "p",
				tip_labels_align = True, 
				node_sizes = 5, 
				scalebar = True )
			toyplot.svg.render( canvas, out_file )

		num_cores = multiprocessing.cpu_count()		
		dir_input = "%s/input_gff" %dir_pan
		dir_temp = "%s/temp_out" %dir_pan
		os.system( "mkdir -p %s" %dir_input )
		file_muscle = "%s/muscle" % path_program
		for genome in os.listdir( dir_genome ) :
			from_gff = "%s/%s/%s.gff" %( dir_genome, genome, genome )
			to_gff = "%s/%s.gff" %( dir_input, genome )
			shutil.copy( from_gff, to_gff )
		
		if len( os.listdir( dir_input ) ) > 1 :
			backbone = "%sphylo_backbone" % path_data
			gff = dir_input + "/*.gff"
			dir_temp0 = "%s/temp_out0" %dir_pan
			dir_temp1 = "%s/temp_out1" %dir_pan
			os.system( "mkdir -p %s" %dir_temp0 )
			os.system( "mkdir -p %s" %dir_temp1 )
			Rtab_file0 = "%s/gene_presence_absence.Rtab" %dir_temp0
			Rtab_file1 = "%s/gene_presence_absence.Rtab" %dir_temp1
			if os.path.exists( Rtab_file0 ) == True and os.path.getsize( Rtab_file0 ) >= 10 * 1024 :
				if os.path.exists( Rtab_file1 ) == True and os.path.getsize( Rtab_file1 ) >= 10 * 1024 :
					pass
				else : 
					panaroo_merged( backbone, dir_temp0, dir_temp1, num_cores )				
			else : 
				panaroo_default( gff, dir_temp0, num_cores )
				panaroo_merged( backbone, dir_temp0, dir_temp1, num_cores )
			out_file = "%s/pan_matrix_input.fasta" %dir_temp1
			out_aln = out_file.replace( ".fasta", ".aln" )
			out_nwk = out_file.replace( ".fasta", ".nwk" )
			df = pd.read_csv( Rtab_file1, sep = "\t" )
			df_transposed = transform_matrix( df )
			create_matrix( df_transposed, out_file )
			if os.path.exists( out_aln ) == True and os.path.getsize( out_aln  ) >= 100 * 1024 : 
				pass
			else : 
				msa_mafft( out_file, out_aln, num_cores )
				make_phylo( file_muscle, out_aln, out_nwk )
			result_phylo = "%s/02.pan_phylo.svg" %dir_out
			generate_tree( out_nwk, result_phylo )
		print ( "PIPeco: Pan phylo process Done..." )