import re
import sys
import vcf

fields = ['Chromosome', 'Position', 'Ref_Allele_Forward_Strand', 'Alt_Allele_Fwd_Strand', 'Average_Depth', 'No_calls', 'Homozygous_calls', 'Het_calls', 'Mutation_Type', 'Mutation_Strength', 'Mutation_Type2', 'Codon_Subst', 'AminoAcid_Subst', 'Locus_Tag', ]
vcf_reader = vcf.Reader(sys.stdin)

first = True

for record in vcf_reader:
	if first:
		for sample in record.samples:
			#for field in ('_gt', '_freq', '_pval'):
			fields.append(sample.sample)

		fields.append("Max_Freq")
		print "\t".join(fields)

		first = False

	result = []

	result.extend([record.CHROM, record.POS, record.REF, record.ALT[0]])

	result.extend([record.INFO['ADP'], record.INFO['NC'], record.INFO['HOM'], record.INFO['HET']])

	parsed_effects = []
	try:
		effects = record.INFO['EFF'].split(",")
		for eff in effects:
			m = re.match("(.*)\((.*)\)", eff)
			if m:
				flags = m.group(2).split("|")
				if 'gene' in flags[7]:
					parsed_effects.extend([m.group(1), flags[0], flags[1], flags[2], flags[3], flags[4]])
	except KeyError, e: pass

	if parsed_effects:
		result.extend(parsed_effects)
	else:
		result.extend([""] * 6)


	freq = 0
	for sample in record.samples:
		if sample.gt_bases:
			bases = sample.gt_bases.split("/")
			if bases[0] == bases[1]:
				if record.REF == bases[0]:
					result.append('')
				else:
					result.append(bases[0])
			else:
				result.append(sample.gt_bases)	
		else:
			result.append('.')
		##result.extend([sample.gt_bases])
		##, sample.data.FREQ, sample.data.PVAL])
		freq = max(freq, sample.data.FREQ)

	result.append(freq)
	
	"""
		maxvars = []
		for col in cols[9:]:
			genotype = col.split(":")

			if col.startswith('1/1'):
				result.append(cols[4])
			elif col.startswith('./.'):
				result.append(".")
			else:
				result.append("%s/%s" % (cols[3], cols[4]))

			try:
				result.append(genotype[6])
				maxvars.append(int(genotype[6].replace("%","")))
			except:
				result.append("")

		if maxvars:
			result.append(str(max(maxvars)))
		else:
			result.append("")

"""
	print "\t".join([str(x) for x in result])


