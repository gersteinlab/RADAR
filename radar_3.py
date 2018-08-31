"""
RADAR version that looks up main scores from precomputed tables
"""

import argparse
import ntpath
import os
import pybedtools

from multiprocessing import Process, Manager, Array, Pool

# files and directories we need
RESOURCES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources/")
MAIN_SCORES_DIR = os.path.join(RESOURCES_DIR, "main_scores/{}/")
ALL_RBP_SITES = os.path.join(RESOURCES_DIR, "all_RBP_peaks_unmerged_labeled_sorted.bed")
SIGNIFICANT_PEAKS = os.path.join(RESOURCES_DIR, "significant_peaks")
MUTATIONAL_BURDEN_MAT = os.path.join(RESOURCES_DIR, "rbp_peak_significance")
REG_POWER_MAT = os.path.join(RESOURCES_DIR, "regulator_pval.txt")

CANCERS = ['BLCA', 'BRCA', 'CESC', 'COAD', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'SKCM', 'PAAD', 'PRAD', 'STAD', 'THCA', 'UCEC']
ASSEMBLIES = ['hg19', 'hg38']

# class that represents a variant and stores its scores
class Variant:
	def __init__(self, key):
		self.key = key
		self.gene_target, self.reg_power, self.mut_burden, self.main = None, None, None, None
	def score_string(self):  # returns tab-separated string of scores in correct order, ends with newline
		tissue_specific = [score for score in (self.gene_target, self.reg_power, self.mut_burden) if score != None]
		if self.main == None:  # if not in regulome
			return "{}\n".format("\t".join(list(self.key) + ["0"] * (6 + len(tissue_specific))))
		return "{}\n".format("\t".join(list(self.key) + self.main + tissue_specific))

def score_string(key, var):  # returns tab-separated string of scores in correct order, ends with newline
	tissue_specific = [score for score in var[1] if score != None]
	if var[0] == None:  # if not in regulome
		return "{}\n".format("\t".join(list(key) + ["0"] * (9 + len(tissue_specific))))
	main_total, ts_total = float(var[0][-1]), sum(int(val) for val in tissue_specific)
	return "{}\n".format("\t".join(list(key) + var[0] + tissue_specific + [str(ts_total), str(main_total + ts_total)]))

# helper function to search main score file for relevant variants
def search_score_files(tup):
	ch, cnt = tup
	curr_found = 0
	variants = dict()
	try:
		with open(os.path.join(MAIN_SCORES_DIR, "{}_scored".format(ch))) as file:
			for line in file:
				line = line.split()
				key = (line[0], line[1], line[2], line[3].upper(), line[4])
				# load main scores for each variant we're considering
				if key in var_set:
					variants[key] = variants[key] = line[5:12]
					curr_found += 1
					if curr_found == cnt:  # break early if we've found all requested variants on this chromosome
						break
	except:
		pass
	return ch, variants

parser = argparse.ArgumentParser(description='RADAR')

parser.add_argument('-b', '--bed', help="Variant BED file", required=True)
parser.add_argument('-a', '--assembly', help="Genome assembly", default="hg19")
parser.add_argument('-o', '--outdir', help="Output directory", required=True)
parser.add_argument('-c', '--cancer', help="Cancer type for tissue-specific scoring", required=False)
parser.add_argument('-kg', '--keygenes', action="store_true", default=False, help="Compute key genes score")
parser.add_argument('-mr', '--mutrec', action="store_true", default=False, help="Compute mutation recurrence score")
parser.add_argument('-rp', '--regpower', action="store_true", default=False, help="Compute RBP regulation power score")

args = parser.parse_args()

assembly = args.assembly
if assembly not in ASSEMBLIES:
	print("Invalid assembly")
	exit(1)
MAIN_SCORES_DIR = MAIN_SCORES_DIR.format(assembly)

# cancer type used for tissue-specific scores
if assembly == 'hg19' and (args.keygenes or args.mutrec or args.regpower):
	cancer_type = args.cancer
	try:
		cancer_index = CANCERS.index(cancer_type)
	except:
		print("Invalid cancer type provided")
		exit(1)

# first preprocess for any tissue-specific scores, which require intersections
variant_string_list = []
with open(args.bed) as file:
	for line in file:
		variant_string_list.append("\t".join(line.split()[:5]))
variants_bedtool = pybedtools.BedTool("\n".join(variant_string_list), from_string=True)

if args.keygenes and assembly == 'hg19':
	# find locations for each requested significant gene in this particular cancer
	significant_peaks_string_list = []
	with open(SIGNIFICANT_PEAKS) as file:
		next(file)  # skip header
		for line in file:
			line = line.split()
			sig = line[4 + cancer_index]
			if sig == "1":  # we only care if this is a significant peak
				significant_peaks_string_list.append("\t".join(line[:3]))
	# intersect variants with significant peaks
	significant_peaks_bedtool = pybedtools.BedTool("\n".join(significant_peaks_string_list), from_string=True)
	variants_bedtool = variants_bedtool.intersect(significant_peaks_bedtool, c=True)

# include mutational burden
if args.mutrec and assembly == 'hg19':
	# string list of burdened RBP peaks
	burdened_peaks_string_list = []
	with open(MUTATIONAL_BURDEN_MAT) as file:
		for line in file:
			line = line.split()
			peak, burdens = line[:3], line[4:]
			if burdens[cancer_index] == "1":  # burdened
				burdened_peaks_string_list.append("\t".join(peak))
	# check how many burdened RBP peaks intersect with each variant
	burdened_peaks_bedtool = pybedtools.BedTool("\n".join(burdened_peaks_string_list), from_string=True)
	variants_bedtool = variants_bedtool.intersect(burdened_peaks_bedtool, c=True)

# include regulatory power
if args.regpower and assembly == 'hg19':
	# load regulatory power matrix
	with open(REG_POWER_MAT, 'r') as file:
		next(file)
		reg_power = set()  # set of RBPs with high regulatory power in requested cancer
		for line in file:
			line = line.split()
			rbp, sig = line[0], line[1 + cancer_index]
			if sig == "1":
				reg_power.add(rbp)
	# check which RBPs each variant intersects
	rbps_bedtool = pybedtools.BedTool(ALL_RBP_SITES)
	variants_bedtool = variants_bedtool.intersect(rbps_bedtool, loj=True)

# list of variants we need to score
variant_list = []  # holds a list of variant keys, same order as given
chromosomes = dict()  # maps chromosomes to number of input variants on that chromosome
variants_ts_scores = dict()

# load in variants to score
for line in variants_bedtool:
	line = list(line)
	key = tuple(line[:5])  # ch, start, stop, ref, alt
	# add any tissue-specific scores requested
	if key not in variants_ts_scores:
		variant_list.append(key)
		chromosomes[key[0]] = chromosomes[key[0]] + 1 if key[0] in chromosomes else 1
		# variants[key] = Variant(key)
		variants_ts_scores[key] = [None, None, None]
		# only need to load these once
		if args.keygenes and assembly == 'hg19':  # gene-target scores are in the 6th column
			val = "1" if line[5] != "0" else "0"  # score is 1 if intersects at least one gene
			variants_ts_scores[key][0] = val
		if args.mutrec and assembly == 'hg19':
			val = "1" if line[5 + (args.keygenes)] != "0" else "0"  # score is 0 if intersects at least one burdened peak
			variants_ts_scores[key][1] = val
	if args.regpower and assembly == 'hg19':
		# rbp in 9th, 10th, or 11th column depending on whether gene-target and mutational burden scores requested
		rbp = line[8 + args.keygenes + args.mutrec]
		# score 1 if intersects at least one powerful rbp in this cancer, else 0
		val = "1" if rbp in reg_power else "0"
		variants_ts_scores[key][2] = val

# search relevant main score files in parallel
var_set = set(variant_list)  # make set for easy lookup
pool = Pool()
inputs = [(ch, chromosomes[ch]) for ch in chromosomes]
variants = dict(pool.map(search_score_files, inputs))

# write scores to output file
head, tail = ntpath.split(args.bed)
var_file = tail or ntpath.basename(head)
bed_name = var_file[:var_file.index(".bed")] if len(var_file) > 4 and var_file[-4:] == ".bed" else var_file
output_file = os.path.join(args.outdir, "{}.radar_out.bed".format(bed_name))
header = ['chr', 'start', 'stop', 'ref', 'alt', 'cross_species_conservation', 'RBP_binding_hub', 'GERP', 'Evofold', 'motif_disruption', 'RBP_gene_association', 'total_universal']
with open(output_file, 'w') as outfile:
	outfile.write("\t".join(header))
	if args.keygenes and assembly == 'hg19': outfile.write("\tkey_genes")
	if args.mutrec and assembly == 'hg19': outfile.write("\tmutation_recurrence")
	if args.regpower and assembly == 'hg19': outfile.write("\tRBP_regulation_power")
	outfile.write("\ttotal_tissue_specific")
	outfile.write("\ttotal_score\n")
	for var in variant_list:
		ch = var[0]
		outfile.write(score_string(var, [variants[ch][var] if var in variants[ch] else None, variants_ts_scores[var]]))
