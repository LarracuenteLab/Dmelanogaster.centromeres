Xiaolu Wei 
Xiaolu_Wei@URMC.Rochester.edu
For figures in Chang et al. 2019 PLoS Biology

TE analysis and distribution along chromosomes
1. grep repeats variants from genome based on annotations in gff file ( “grep_repeat_gff_genome.py” )
   Alternatively, we can also blast repeat consensus sequence to genome and parse the output to get repeat variants. ("parse_blast_output_consensus_genome.py")
2. blast variants sequences to consensus sequence.
3. Parse blast output. (“parse_blast_output_variant_consensus.py”). It generates _blast.summary and _plot.summary
   Note: when parse blast output, try the two different ways to interpret the aligment in query (consensus) in the script, see if they agrees in deciding if the repeat is full length.
   Note: To do phylogeny analysis, output big elements (>100bp) and align them in Geneious to get alignment file. Then run Raxml to get tree file, and parse it using "parse_raxml_output.py" to generate input for visualization in R ("plot_raxml_tree.R").
4. summary all repeats by contigs and chromosomes. (“summarize_repeat_chromosome_position.sh" )
5. summary all repeats on chromosomes by contig order on chromosomes. (“summarize_repeat_chromosome_position.py")
6. plot in R. ("Fig4_FigS11.3.5.19.R")