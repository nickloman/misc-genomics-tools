from Bio import SeqIO
import sys

fasta = open(sys.argv[1])
headers = [arg for arg in sys.argv[2:]]

for rec in SeqIO.parse(fasta, "fasta"):
	if rec.id in headers:
		SeqIO.write([rec], sys.stdout, "fasta")	
