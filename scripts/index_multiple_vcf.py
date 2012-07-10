
#import pysam
#import sys

# usage: <all variant positions.vcf> <each individual sample>

#all_variants_vcf_fn = sys.argv[1]

#vcf = pysam.VCF()
#for ln in vcf.parse(open(all_variants_vcf_fn)):
#	print ln

#for ln in vcf.fetch():
#	print ln

import vcf
import sys
import subprocess
import pickle
from tabix import Tabix

class Output:
	def __init__(self, samples, results):
		self.samples = samples
		self.results = results


positions = set()

for ln in open(sys.argv[1]):
	fn = ln.rstrip()
	vcf_reader = vcf.VCFReader(open(fn), 'rb')
	for record in vcf_reader:
		positions.add(( record.CHROM, record.POS ))


samples = []
results = {}

for ln in open(sys.argv[2]):
	sample, fn = ln.rstrip().split("\t")
	samples.append(sample)
	results[sample] = {}

#	p = subprocess.Popen(['bcftools', 'view', fn], stdout=subprocess.PIPE)
#	vcf_reader = vcf.VCFReader(p.stdout, 'rb')
#	for record in vcf_reader:
#		pos = (record.CHROM, record.POS)

	print fn
	tab = Tabix(fn)
	for pos in positions:
		search = "%s:%d-%d;" % (pos[0], pos[1], pos[1])
		try:
			itr = tab.fetch(search)
			rec = itr.next()
		except StopIteration:
			print "can't find: %s" % (search,)
			pass

		cols = rec.split("\t")

		record = {}
		record['REF'] = cols[3]
		record['ALT'] = cols[4]
		record['QUAL'] = float(cols[5])
		flags = cols[7].split(";")
		for f in flags:
			try:
				key, val = f.split("=")
				if key == 'AF1':
					record['AF1'] = float(val)
				if key == 'MQ':
					record['MQ'] = int(val)
				if key == 'DP':
					record['DP'] = int(val)
			except:
				pass

		# look it up
		results[sample][pos] = record

pickle.dump(results, open(sys.argv[3], "w"))
