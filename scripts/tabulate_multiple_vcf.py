
import pickle
import sys

def get_base( r ):
	if r['ALT'] != '.':
		return r['ALT']
	else:
		return r['REF']

results = pickle.load(open(sys.argv[1]))

samples = results.keys()
n_samples = len(samples)
positions = set()
sample_bases = {}
for s in samples:
	for tup in results[s].keys():
		positions.add(tup)
	sample_bases[s] = []

reference_bases = []

for n, pos in enumerate(sorted(positions)):
	print n, pos
	results_array = []
	for s in samples:
		if pos in results[s]:
			results_array.append ( results[s][pos] )

	print results_array

	# loosk like an indel
	if [ r for r in results_array if len(r['REF']) > 1 ]:
		continue

	if [ r for r in results_array if len(r['ALT']) > 1 ]:
		continue

	if len(results_array) != n_samples:
		continue

	# mapping quality - could do this per base
	if [ r for r in results_array if 'MQ' in r and  r['MQ'] < 30 ]:
		continue

	# homozygote

	if [ r for r in results_array if 'AF1' in r and r['AF1'] > 0.0 and r['AF1'] < 1.0]:
		continue

	reference_bases.append(results_array[0]['REF'])

	for s in samples:
		sample_bases[s].append( get_base( results[s][pos] ) )

print ">ref"
print "".join(reference_bases)

for s in samples:
	print ">%s" % (s,)
	print "".join(sample_bases[s])
