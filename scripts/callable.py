import vcf
import sys
import numpy


vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))
itr = vcf_reader.__iter__()
record = itr.next()

labels = [sample.sample for sample in record.samples]

samples = numpy.array([ 0 for sample in record.samples ])

def one_if_not_null(x):
	if x is None:
		return 0
	else:
		return 1

positions = 0

while True:
	samples += [one_if_not_null(sample.gt_bases) for sample in record.samples]
	positions += 1

	try:
		record = itr.next()
	except StopIteration:
		break

print "#positions: %d" % (positions,)
for n, label in enumerate(labels):
	print "%s,%s,%.02f" % (label, samples[n], 100.0 / positions * samples[n])

