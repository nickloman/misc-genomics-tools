import vcf
import sys
import re

r = re.compile(sys.argv[2])
vcf_reader = vcf.Reader(filename=sys.argv[1])
vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
for record in vcf_reader:
	if any(r.match(str(s)) for s in record.INFO.values()):
        	vcf_writer.write_record(record)
