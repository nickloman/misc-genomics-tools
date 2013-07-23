import cPickle
from Bio import SeqIO
import sys
from SeqUtils import get_snp_mutation_effect, OneBasedSystem, get_seq_1_based

sequence_maps = cPickle.load(open(sys.argv[1]))

def get_overlapping_features(rec, pos):
    for feature in rec.features:
        if feature.type == 'source':
            continue
        if feature.type not in ['CDS', 'tRNA', 'rRNA']:
            continue
        if feature.location.nofuzzy_start+1 <= pos and feature.location.nofuzzy_end >= pos:
            yield feature

class FragmentedGenome:
    def __init__(self, recs):
        self.recs = list(recs)

    def translate_contig_coords_to_concatenated_coords(self, contig, pos):
        n = 0
        for r in self.recs:
            if r.id == contig:
                return n + int(pos)
            n += len(r)
        print contig, pos

def distances_to_alignment_break(sequencemap, pos):
    fwd = 0
    rev = 0
    for n, c in enumerate(sequencemap[pos:]):
        if c == '-' or c == ' ':
            fwd = n
            break
    revstr = sequencemap[:pos]
    revstr = revstr[::-1]
    for n, c in enumerate(revstr):
        if c == '-' or c == ' ':
            rev = n
            break
    return fwd, rev
        

def go():
    g = FragmentedGenome(SeqIO.parse(open(sys.argv[2]), "fasta"))
    reference_rec = SeqIO.parse(open(sys.argv[3]), "genbank").next()

    for n, ln in enumerate(sys.stdin):
        if not n: continue

        cols = ln.rstrip().split("\t")

        pos = g.translate_contig_coords_to_concatenated_coords(cols[0], cols[1])

        fwd, rev = distances_to_alignment_break(sequence_maps['seq2_base_lookup'], pos-1)

        reference_base = sequence_maps['seq2_coord_lookup'][pos-1]
        if type(reference_base) is int:
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("\t".join(cols), pos, sequence_maps['seq2_base_lookup'][pos-1], fwd, rev, reference_base + 1, sequence_maps['seq2_strand_lookup'][pos-1]),
        else:
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("\t".join(cols), pos, '-', '-', '-', '-', '-'),

        if type(reference_base) is not int:
            print
            continue

        seenset = set()
        for snp in cols[2:]:
            if sequence_maps['seq2_base_lookup'][pos-1].lower() == snp.lower():
                continue

            if sequence_maps['seq2_strand_lookup'][pos-1] == '-':
                strand_snp = dict(zip('ACGTacgt','TGCAtgca'))[snp]
            else:
                strand_snp = snp

            cds = [f for f in get_overlapping_features(reference_rec, reference_base + 1) if f.type == 'CDS']
            if cds:
                cds = cds[0]
                snp_type, \
                amino_acid_reference, \
                amino_acid_consensus, \
                substitution_score, \
                ref_base, \
                ref_sequence, \
                cons_sequence = \
                    get_snp_mutation_effect(str(reference_rec.seq),
                        OneBasedSystem(reference_base + 1),
                        strand_snp,
                        OneBasedSystem(cds.location.nofuzzy_start + 1),
                        OneBasedSystem(cds.location.nofuzzy_end),
                        cds.strand)
                if strand_snp.lower() != ref_base.lower() and snp not in seenset:
#                    print "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (cds.qualifiers['locus_tag'][0], cds.qualifiers['product'][0], ref_base, snp_type, amino_acid_reference, amino_acid_consensus, substitution_score, ref_sequence, cons_sequence),
                    print "\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (cds.qualifiers['locus_tag'][0], cds.qualifiers['product'][0], ref_base, snp_type, amino_acid_reference, amino_acid_consensus, substitution_score),
                    seenset.add(snp)
        print    

go()
