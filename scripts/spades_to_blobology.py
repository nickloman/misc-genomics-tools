# usage:
# python spades_to_blobology.py assembly1 assembly2 assembly3 > stats.txt
# R script
# st<-read.table("stats.txt", header=TRUE)
# ggplot(st, aes(y=Coverage, x=GC, size=Length)) + geom_point() + facet_wrap( ~ Sample )

from Bio import SeqIO
from Bio.SeqUtils import GC
import sys

def contig_stat(fn, rec):
	cols = rec.id.split("_")
	return "%s\t%s\t%s\t%s\t%s" % (fn, rec.id, len(rec), GC(rec.seq), cols[5])

def process(fn):
	print "\n".join([contig_stat(fn, rec) for rec in SeqIO.parse(open(fn), "fasta")])

print "Sample\tContig\tLength\tGC\tCoverage"
for fn in sys.argv[1:]:
	process (fn)
