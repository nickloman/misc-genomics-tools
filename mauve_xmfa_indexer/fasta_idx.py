
from Bio import SeqIO
import sys
import pickle
import re

sequencemaps = {
    'seq1_base_lookup' : list(' ' * 10000000),
    'seq2_base_lookup' : list(' ' * 10000000),
    'seq1_coord_lookup' : list(' ' * 10000000),
    'seq2_coord_lookup' : list(' ' * 10000000),
    'seq1_strand_lookup' : list(' ' * 10000000),
    'seq2_strand_lookup' : list(' ' * 10000000)
}

def parse_hdr(rec):
    hdr = {}
    m = re.search('(\d+):(\d+)-(\d+) ([+-]) (.*)', rec.description)
    hdr['seqnum'] = int(m.group(1))
    hdr['start'] = int(m.group(2))
    hdr['end'] = int(m.group(3))
    hdr['strand'] = m.group(4)
    return hdr

def crunch(mapid, rec1, rec2):
    print >>sys.stderr, len(rec1), len(rec2)

    h1 = parse_hdr(rec1)
    h2 = parse_hdr(rec2)

    # if h1['strand'] != '+' or h2['strand'] != '+':
    #    return

    if h1['strand'] == '+':
        seq1_absolute_position = h1['start'] - 1
    else:
        seq1_absolute_position = h1['end'] - 1

    if h2['strand'] == '+':
        seq2_absolute_position = h2['start'] - 1
    else:
        seq2_absolute_position = h2['end'] - 1

    for seq1_alignment_position, c in enumerate(rec1.seq):
        if c != '-':
            # print >>sys.stderr, c, rec2.seq[n]

            sequencemaps["%s_base_lookup" % (mapid,)][seq1_absolute_position] = str(rec2.seq[seq1_alignment_position])
            sequencemaps["%s_strand_lookup" % (mapid,)][seq1_absolute_position] = h2['strand']

            if rec2.seq[seq1_alignment_position] != '-':
                sequencemaps["%s_coord_lookup" % (mapid,)][seq1_absolute_position] = seq2_absolute_position

            if h1['strand'] == '+':
                seq1_absolute_position += 1
            else:
                seq1_absolute_position -= 1

        if rec2.seq[seq1_alignment_position] != '-':
            if h2['strand'] == '+':
                seq2_absolute_position += 1
            else:
                seq2_absolute_position -= 1

def go():
    aln = []
    for rec in SeqIO.parse(open(sys.argv[1]), "fasta"):
        print >>sys.stderr, rec.description

# i> 1:1-47968 - data/reference_sequences/PAO1_NC_002516.gbk

        if "=" == str(rec.seq)[-1]:
            end_of_alignment = True
            rec.seq = rec.seq[0:-1]
        else:
            end_of_alignment = False

        aln.append(rec)
        if end_of_alignment:
            if len(aln) == 2:
                crunch('seq1', aln[0], aln[1])
                crunch('seq2', aln[1], aln[0])
            aln = []

    pickle.dump(sequencemaps, sys.stdout)

go()

