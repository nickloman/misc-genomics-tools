import vcf
import sys
import re
from collections import defaultdict
from optparse import OptionParser
import operator

def read_table(fn):
	table = []
	colnames = []
	headers = False
	for ln in open(fn):
		ln = ln.rstrip()
		if ln.startswith('#') and headers == False:
			colnames = ln[1:].split("\t")
			headers = True
			continue
		if ln.startswith('#'):
			continue
		cols = ln.split("\t")
		d = {}
		for key in colnames:
			d[key] = ''
			for n, v in enumerate(cols):
				try:
					d[colnames[n]] = v
				except:
					print >>sys.stderr, "Warning misformed %s" % (ln,)
		table.append(d)
	return colnames, table

def go(args, options):
	metadata_columns, metadata = read_table(options.metadata)
	fh_alleles = open("alleles.txt", "w")
	fh_samples = open("samples.txt", "w")

	alleles = defaultdict(str)
	nocalls = defaultdict(int)
	nocalls_by_rec = defaultdict(list)
	nocalls_list = []

	vcf_reader = vcf.Reader(open(args[0]), "r")
	for record in vcf_reader:
		for sample in record.samples:
			if sample.gt_bases is None:#
				nocalls[sample.sample] += 1

				nocalls_by_rec["%s-%s" % (record.CHROM, record.POS)].append(sample)

				nocalls_list.append([record.CHROM, record.POS, record.REF, record.ALT, sample.sample, sample.data.SDP])

		bases = [sample.gt_bases for sample in record.samples if sample.gt_bases is not None]
	
		if len(bases) != len(record.samples):
			continue

		incorrect_lengths = [sample for sample in record.samples if len(sample.gt_bases.split("/")[0]) != 1]

		if incorrect_lengths: ## e.g. indels
			continue	

		for sample in record.samples:
			alleles[sample.sample] += sample.gt_bases.split("/")[0]
	
	unique_alleles = set(alleles.values())
	allele_lookup = dict([(allele, n+1) for n, allele in enumerate(unique_alleles)])

	print metadata

	print >>fh_samples, "Sample\tProfile\t",
	print >>fh_samples, "\t".join(metadata_columns)
	for sample, genotype in alleles.iteritems():
		print >>fh_samples, "%s\t%s\t" % (sample, allele_lookup[genotype]) ,
		print sample
		for sample_metadata in metadata:
			print "  " + sample_metadata['Sample']
			if sample_metadata['Sample'] == sample:
				print >>fh_samples, "\t".join([sample_metadata[key] for key in metadata_columns if key != 'Sample']) ,
		print >>fh_samples

	for allele, allele_number in allele_lookup.iteritems():
		print >>fh_alleles, ">%s\n%s" % (allele_number, allele)

	for sample, value in sorted(nocalls.iteritems(), key=operator.itemgetter(1)):
		print "%s => %s" % ( sample, value )

	for loc, samplelist in sorted(nocalls_by_rec.iteritems(), key=operator.itemgetter(1)):
		print "%s: " % (loc,)
 		for sample in samplelist:
			print "   %s" % (sample)

	if options.nocallsfile:
		nocall_fh = open(options.nocallsfile, "w")
		print "CHROM\tPOS\tREF\tALT\tSample\tSDP"
		for nocall in nocalls_list:
			print >>nocall_fh, "\t".join([str(s) for s in nocall])
		nocall_fh.close()

def main():
	usage = "usage: %prog [options] vcffile"
	parser = OptionParser(usage)
	parser.add_option("-m", "--metadata", dest="metadata",
					  help="load metadata from METADATA")
	parser.add_option("-s", "--samples", dest="samplefile",
					  help="output sample info to SAMPLEFILE (default: <vcffile>_samples.txt")
	parser.add_option("-a", "--alleles", dest="allelesfile",
					  help="output alleles file to ALLELESFILE (default: <vcffile>_alleles.txt")
        parser.add_option("-n", "--nocalls", dest="nocallsfile",
                                          help="output alleles file to NOCALLSFILE (default: <vcffile>_nocalls.txt")
	parser.add_option("-v", "--verbose",
					  action="store_true", dest="verbose")
	parser.add_option("-q", "--quiet",
					  action="store_false", dest="verbose")

	(options, args) = parser.parse_args()
	if len(args) != 1:
		parser.error("incorrect number of arguments")
	if options.verbose:
		print "reading %s..." % options.filename

	go(args, options)

if __name__ == "__main__":
	main()

