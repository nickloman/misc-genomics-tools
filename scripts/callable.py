import vcf
import sys
import numpy

def one_if_not_null(x):
        if x is None:
                return 0
        else:
                return 1

vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))
itr = vcf_reader.__iter__()
record = itr.next()

labels = [sample.sample for sample in record.samples]
samples = numpy.array([ 0 for sample in record.samples ])

calls = nocalls = 0

while True:
	mx = [one_if_not_null(sample.gt_bases) for sample in record.samples]
	if 0 in mx:
		nocalls += 1
	else:
		calls += 1

	samples += mx

	try:
		record = itr.next()
	except StopIteration:
		break

print "#calls: %d, nocalls: %d" % (calls, nocalls,)
print "Label,Calls,Perc"
for n, label in enumerate(labels):
	print "%s,%s,%.02f" % (label, samples[n], 100.0 / (calls+nocalls) * samples[n])

