import sys
from Bio import SeqIO

recs = list(SeqIO.parse(open(sys.argv[1]), "fasta"))

ref = recs[0]
rest = recs[1:]

snp = 1
print "Sample\tGenomePos\tSNPNum\tVariant"
for n in xrange(len(ref)):
	var = 0
	if len([sample for sample in rest if sample[n] != ref[n]]) == len(rest):
		continue

	for sample in rest:
		if sample[n] == ref[n]:
			print "%s\t%s\t%s\t%s" % (sample.id, n+1, snp, 0)
		else:
			print "%s\t%s\t%s\t%s" % (sample.id, n+1, snp, 1)

	snp += 1

