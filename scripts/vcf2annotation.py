import re
import sys
import vcf
from collections import defaultdict
from renamefasta import getlabels

##Label    Effect  Effect
##Sample   1       0
##Sample   1       0 

vcf_reader = vcf.Reader(sys.stdin)
samples = defaultdict(list)
effects = []

for record in vcf_reader:
	effect = record.INFO['EFF'].split(",")
	for eff in effect:
		m = re.match("(.*)\((.*)\)", eff)
		if m:
			flags = m.group(2).split("|")
			effects.append("%s_%s" % (flags[4], flags[3]))

	for sample in record.samples:
		samplename = sample.sample
		if sample.gt_bases:
			bases = sample.gt_bases.split("/")
			if bases[0] == bases[1]:
				if record.REF == bases[0]:
					samples[samplename].append('1')	
				else:
					samples[samplename].append('2')
			else:
				samples[samplename].append('3')
		else:
			samples[samplename].append('0')

labels = getlabels(sys.argv[1])

print "taxa\t" + "\t".join(effects)
for s, genotypes in samples.iteritems():
	print labels[s] + "\t" + "\t".join(genotypes)
