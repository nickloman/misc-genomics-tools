import subprocess
import vcf
import random
import sys

from reportlab.platypus import BaseDocTemplate, SimpleDocTemplate, Paragraph, Spacer, Preformatted, Frame, PageTemplate, KeepTogether
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.pagesizes import letter, portrait, landscape, A4
from reportlab.lib.units import inch

PAGE_HEIGHT=defaultPageSize[0]; PAGE_WIDTH=defaultPageSize[1]
styles = getSampleStyleSheet()

def get_vcf(stream):
	return [(v.CHROM, v.POS) for v in vcf.Reader(stream)]

monkeys = random.sample(get_vcf(sys.stdin), 18)

def get_ttview(chrom, pos, rows):
	X = 100
	args = ["../samtools-utilities/bin/ttview",
		"-X", str(X),
		"-g", "%s:%s" % (chrom, max(0, pos - (X/2))),
		"../demodata/sample_906_sorted.bam",
		"../demodata/Pseudomonas_aeruginosa_NCGM2S1.fasta"]
	p = subprocess.Popen(args, stdout=subprocess.PIPE)
	out, err = p.communicate()
	return "\n".join(out.split("\n")[0:rows])

Title = "Hello world"
pageinfo = "platypus example"
def myFirstPage(canvas, doc):
	pass

def go():
	#create the basic page and frames
	doc = BaseDocTemplate("phello.pdf", leftMargin = 10, rightMargin = 0, topMargin = 0, bottomMargin = 0)
	doc.pagesize = landscape(A4)

	frameCount = 2
	frameWidth = doc.height/frameCount
	frameHeight = doc.width-.05*inch
	frames = []
	#construct a frame for each column
	for frame in range(frameCount):
		leftMargin = doc.leftMargin + frame*frameWidth
		column = Frame(leftMargin, doc.bottomMargin, frameWidth, frameHeight)
		print leftMargin, doc.bottomMargin, frameWidth, frameHeight

		frames.append(column)
	template = PageTemplate(frames=frames)
	doc.addPageTemplates(template)

#	doc = SimpleDocTemplate("phello.pdf", id='TwoColumns')
#	doc.pagesize = landscape(A4) ## TODO: make configurable

	Story = []
	style = styles["Normal"]
	style.fontName = 'Courier'
	style.fontSize = 6
	for monkey in monkeys:
		p = Preformatted(get_ttview(monkey[0], monkey[1], 15), style)
		Story.append(KeepTogether(p))
		Story.append(Spacer(1,0.05*inch))
	doc.build(Story)

#, onFirstPage=myFirstPage)

go()
