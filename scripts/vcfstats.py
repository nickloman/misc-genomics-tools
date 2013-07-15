### Ensure your SNP caller will produce het calls,e.g.
### java -jar bin/VarScan.jar mpileup2snp --min-coverage 5 --min-var-freq 0.1 --p-value 0.005  --output-vcf

import vcf
import sys
from collections import Counter, defaultdict

counts = defaultdict(Counter)

vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))
for record in vcf_reader:
	for sample in record.samples:
		counts[sample.sample].update([sample['GT']])

print "Sample\tHET\tHOM"
for sample, counts in counts.iteritems():
	print sample,
	for gt in ['0/1', '1/1']:
		print "\t" + str(counts[gt]),
	print 
