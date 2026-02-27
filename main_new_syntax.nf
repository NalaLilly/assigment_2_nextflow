	params {
    inputFile: Path = 'data/bacterial_dna.fasta'
	cutoff: Integer = 60
}



process countGC {

    input:
    path inputFile
    val cutoff
	
	output:
    path 'output.txt'
    
    script:
    """
Rscript -e '
	
	cutoff <- as.integer("${cutoff}")
	fasta_file <-"${inputFile}"
	
	library("Biostrings")
	seqs <- readDNAStringSet(fasta_file)
	sequance <- as.character(seqs)
	GC_count <- letterFrequency(seqs, letters = "GC")
	SeqID <- names(seqs)

	result_file = "output.txt"
	cat("",file = result_file)
						 
	for(i in 1:length(sequance)) {
	  if(GC_count[i] > (cutoff)) {

		cat(SeqID[i],"\n",sequance[i], "\n", "\n",file = result_file, append = TRUE)
	}}
	'
	"""
}
	
workflow {
	
	main:
	
	//path for the seq fasta file fasta
    seq_file_ch = channel.fromPath(params.inputFile)
	

    // GC count in a sequancer
	countGC(seq_file_ch, params.cutoff)
    
	publish : 
	first_output = countGC.out
}
	
output{
	first_output{
		path 'sequance'
		mode 'copy'
	}
}

