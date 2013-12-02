from Bio import SeqIO
import sys

def getlabels(fn):
	relabel = {}
	unknownn = 1
	for ln in open(fn):
		cols = ln.rstrip().split("\t")
		try:
			cols[1] = cols[1].replace(' ', '_')
			if cols[1] in relabel.values():
				cols[1] += '_'
			relabel[cols[0]] = cols[1]
		except:
			relabel[cols[0]] = 'unknown%d' % (unknownn)
			unknownn += 1
	return relabel

if __name__ == "__main__":
	relabel = getlabels(sys.argv[1])

	for rec in SeqIO.parse(sys.stdin, "fasta"):
		if rec.id in relabel:
			rec.id = relabel[rec.id]
			rec.description = ''
		SeqIO.write([rec], sys.stdout, "fasta")
	

