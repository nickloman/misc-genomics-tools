from csv import DictReader
from Bio.Phylo import NexusIO
from Bio.Phylo.BaseTree import Clade
import sys

def get_comments(keys):
	# ! = boolean
	# & = number
	# string
	return ",".join(["&" + k + "=\"" + v + "\"" for k, v in keys if k != "" and v != ""])

d = {}

csv = DictReader(open(sys.argv[2]), dialect='excel-tab')
for row in csv:
	d[row['SampleID']] = row	

tree = list(NexusIO.parse(open(sys.argv[1])))[0]
for t in tree.find_clades():
	if t.name:
		t_name = t.name.strip("'")
		if t_name in d:
			s = d[t_name]
			keys = [[k, s[k]] for k in sorted(s.keys()) if k != "#SampleID"]
			t.comment = "%s" % (get_comments(keys))

NexusIO.write([tree], sys.stdout)

#begin taxa;
#       CD72
#       CD26[&!name="big monkey"][&!jeff="hello"]
#       CD125

"""
taxa = 0
for ln in open(sys.argv[1]):
	if ln.startswith('begin taxa'):
		taxa = 1	
	if taxa and 'taxlabels' in ln:
		taxa = 2
	if taxa == 2:
		sample = ln.strip()
		if sample in d:
			sys.stdout.write("\t%s" % (sample))
			s = d[sample]
			keys = [k for k in sorted(d.keys()) if k != "SampleID"]
			for key in keys:
				sys.stdout.write("&!%s=\"%s\"" % (key, s[key]))
			sys.stdout.write("\n")
			continue
	sys.stdout.write(ln)
"""
