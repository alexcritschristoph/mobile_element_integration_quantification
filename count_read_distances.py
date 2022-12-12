import glob
import pysam
import pandas as pd
import argparse
import numpy as np
from Bio import SeqIO
from collections import defaultdict


def main(bam_directory, 
	contig_name,
	mobile_element_start, 
	mobile_element_stop, 
	read_length = 150,
	max_insert_length = 500):

	results = defaultdict(list)

	# For each mapping file
	for fn in glob.glob(bam_directory.rstrip("/") + "/*.bam"):
		sample =fn.split("/")[-1].split(".bam")[0] # get sample name
		print(sample)
		samfile = pysam.AlignmentFile(fn, "rb") # open that mapping BAM file
		
		all_reads = 0
		circular_reads = 0 # Circular reads
		virus_to_host_reads = 0 # Reads from integrated element
		host_to_host_reads = 0 # Reads from host genome when element is exised
		insert_sizes = []

		## Get coverage of host
		# get genome's coverage from the start to the START of the integrated element
		genome1 = samfile.count_coverage(contig_name, start=0, stop=mobile_element_start) 
		genome1_cov = np.sum([genome1[0], genome1[1], genome1[2], genome1[3]], axis=0) # sum A, C, G, T

		# get genome's coverage from the STOP of the integrated element to the end of the genome
		genome2 = samfile.count_coverage(contig_name, start=mobile_element_stop)
		genome2_cov = np.sum([genome2[0], genome2[1], genome2[2], genome2[3]], axis=0) # sum A, C, G, T

		# get mean coverage of all sites in the genome (both halves)
		genome_cov = round(np.mean(list(genome1_cov) + list(genome2_cov)), 2)

		# # get coverage of mobile element
		mini1 = samfile.count_coverage(contig_name, start=mobile_element_start, stop=mobile_element_stop)
		mini1_cov = np.mean(np.sum([mini1[0], mini1[1], mini1[2], mini1[3]], axis=0)) # sum A, C, G, T
		mini_cov = round(mini1_cov, 2) # easy (no halves)!

		### Determine mean insert size from first 10 Kb of the contig
		## Insert size here means distance from R1 start to R2 stop - will be 0 if reads completely overlap.
		for pileupcolumn in samfile.pileup(contig_name, 1000, 10000, 
				truncate=True,min_mapping_quality=10, ignore_overlap=False):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					if pileupread.alignment.is_reverse:
						insert_sizes.append(pileupread.alignment.reference_start - pileupread.alignment.next_reference_start)
		mean_insert_size = round(np.nan_to_num(np.mean(insert_sizes)))

		### Find all viral circular and integrated reads
		for pileupcolumn in samfile.pileup(contig_name, mobile_element_start, mobile_element_start+max_insert_length, 
				truncate=True,min_mapping_quality=10, ignore_overlap=False):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					## If the read 
					if pileupread.alignment.reference_start > mobile_element_start and pileupread.alignment.reference_end > mobile_element_start:
						if pileupread.alignment.is_reverse: ## Look at reverse reads going OUT from the start of the sequence							
							all_reads += 1
							if pileupread.alignment.next_reference_start > mobile_element_stop - max_insert_length:
								circular_reads += 1
							elif pileupread.alignment.next_reference_start < mobile_element_start - (read_length / 2.0):
								virus_to_host_reads += 1

		# Find all host to host reads (from genomes where the element has been excised)
		for pileupcolumn in samfile.pileup(contig_name, mobile_element_start - max_insert_length, mobile_element_start, 
				truncate=True,min_mapping_quality=10, ignore_overlap=False):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					## If the read 
					if pileupread.alignment.reference_start < mobile_element_start and pileupread.alignment.reference_end < mobile_element_start:
						if pileupread.alignment.is_forward:
							if pileupread.alignment.next_reference_start > mobile_element_stop:
								host_to_host_reads += 1
		results['Sample'].append(sample)
		results['Genome_coverage'].append(genome_cov)
		results['Mge_coverage'].append(mini_cov)
		results['Integrated_reads'].append(virus_to_host_reads)
		results['Circular_reads'].append(circular_reads)
		results['Excised_reads'].append(host_to_host_reads)

	results = pd.DataFrame(results)
	results.to_csv("mobile_element_reads_quantification.tsv", sep="\t", index=None)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Pipeline for detecting Illumina read pairs that can be used to infer if a prophage or mobile element is active, integrated, or excised from a host genome.')

	parser.add_argument("-b", "--bam_file_directory", action="store", default=None, required=True, \
		help='Directory of cleaned and processed FASTQ files')

	parser.add_argument("-c", "--contig", action="store", default=None, required=True, \
		help='Contig containing the integrated mobile element.')

	parser.add_argument("-s", "--start", action="store", default=None, required=True, \
		help='Starting position (bp) of the integrated mobile element / prophage in the contig.')

	parser.add_argument("-e", "--end", action="store", default=None, required=True, \
		help='Ending position (bp) of the integrated mobile element / prophage in the contig.')

	parser.add_argument("-r", "--read_length", action="store", default='150', \
		help='Sequencing read length (either 150 or 75 bp)')

	parser.add_argument("-i", "--max_insert_length", action="store", default='500', \
		help='max_insert_length')


	args = parser.parse_args()

	main(bam_directory = args.bam_file_directory,
		contig_name = args.contig,
		mobile_element_start = int(args.start),
		mobile_element_stop = int(args.end),
		read_length = int(args.read_length),
		max_insert_length = int(args.max_insert_length))