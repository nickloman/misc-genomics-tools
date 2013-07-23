from Bio.Seq import Seq, reverse_complement, translate
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC

import tempfile, os
import sys

def translate_with_x(seq, remove_stop=False):
    #standard_table = makeTableX(CodonTable.ambiguous_dna_by_id[1])
    #translator = Translator(standard_table)
    #seq = translator.translate(seq)
    seq = seq.translate(table=1)
    if remove_stop and seq.endswith('*'):
        return seq[:-1]
    return seq

def bacterial_translate(seq):
    return seq.translate(table=11)

# This is for stuff out of 
def get_seq_0_based(data, start, end, strand):
    # Repeat after me:
    #   biopython is 0-based
    #   biosql is 1-based
    #   python is 0-based (slices)
    # so here, we don't need to do anything

    start = int(start)
    end = int(end)
    seq = data[start:end]
    if strand:
        if int(strand) == -1:
            # we might have unicode here hence the cast
            seq = reverse_complement(str(seq))
        elif int(strand) == 1:
            pass
        else:
            print "Unknown strand type: %s" % (feature.strand)
    return seq

def get_seq_1_based(data, start, end, strand):
    # Repeat after me:
    #   biosql is 1-based
    #   python is 0-based (slices)
    # so here, we need to interpret stuff

    start = int(start) - 1
    end = int(end)
    seq = data[start:end]
    if strand:
        if int(strand) == -1:
            # we might have unicode here hence the cast
            seq = reverse_complement(str(seq))
        elif int(strand) == 1:
            pass
        else:
            print "Unknown strand type: %s" % (strand)
    return seq

def GC(seq):
    gc = d.get('G',0) + d.get('C',0) + d.get('g',0) + d.get('c',0) + d.get('S',0) + d.get('s',0)

    if gc == 0: return 0
    return gc*100.0/len(seq)

def simple_alignment(seq1, seq2):
    align = ""
    seq2_len = len(seq2)
    for i in xrange(0, len(seq1)):
        if seq2_len <= i:
            align = ''
        elif seq1[i] == seq2[i]:
            align += '*'
        else:
            align += '-'
    return align

def mygetattr(list, feature, default):
    f = list.get(feature, [default])
    return f[0]

def fasta_identifier(qualifier_dict, seqfeature_id=0, genome_id=None):
    gene = qualifier_dict.get('gene', '')
    if gene.find(' ') != -1:
       gene = ''
    
    id = "lcl|"
    if genome_id:
        id += str(genome_id) + "|"
    else:
        id += "|"

    return "%s%s|%s|%s" %  \
        (id,
         qualifier_dict.get('locus_tag', qualifier_dict.get('systematic_id', '')),
         gene,
         seqfeature_id)

def identifier(dict):
    if "locus_tag" in dict:
        return dict['locus_tag']
    if "systematic_id" in dict:
        return dict['systematic_id']
    raise ValueError

def flatten(dict):
    newdict = {}
    for x, y in dict.iteritems():
        newdict[x] = y[0]
    return newdict

def best_identifier(qualifiers):
    return fasta_identifier(flatten(qualifiers))

def get_gene(record, feature, offset = 0):
    gene = {}
    gene['identifier'] = fasta_identifier(flatten(feature.qualifiers))
    gene['description'] = mygetattr(feature.qualifiers, 'product', '')

    if offset:
        start = max(feature.location.nofuzzy_start - offset, 0)
        end = min(feature.location.nofuzzy_end + offset, len(record.seq.data))
        gene['sequence'] = get_seq_0_based(record.seq.data, start, end, feature.strand)
    else:
        gene['sequence'] = get_seq_0_based(record.seq.data, feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.strand)
    return gene

def get_protein(record, feature):
    protein = {}
    protein['identifier'] = fasta_identifier(flatten(feature.qualifiers))
    protein['description'] = mygetattr(feature.qualifiers, 'product', '')

    # prefer the annotated translation versus our own one
    if feature.qualifiers.has_key('translation'):
        # what if more than one translation?
        protein['sequence'] = feature.qualifiers['translation'].pop(0)
    else:
        protein['sequence'] = translate(get_seq_0_based(record.seq.data, feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.strand))
    if not len(protein['sequence']):
        print >>sys.stderr, "could not translate %s" % (protein['identifier'])

    return protein

def translate_proteins(record):
    seqs = []
    for feature in record.features:
        seq = Seq(translate(get_seq_0_based(record.seq.data, feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.strand)), IUPAC.IUPACProtein)
        seqs.append(SeqRecord(seq, id=identifier(flatten(feature.qualifiers)), description=mygetattr(feature.qualifiers, 'product', '')))
    return seqs

def pseudo(feature):
    # check if pseudo
    if feature.qualifiers.has_key('pseudo'):
        return True
    return False

def write_proteins(handle, record):
    proteins = [get_protein(record, feature) for feature in record.features if feature.type == 'CDS' and not pseudo(feature)]
    for protein in proteins:
        handle.write(">%d|%s\n" % (bioentry_id, protein['identifier']))
        for n in xrange(0, len(protein['sequence']), 60):
            handle.write(protein['sequence'][n:n+60])
            handle.write("\n")
        handle.write("\n")

def dump_to_temporary_file(seqrecords, fmt='fasta'):
    fh, fn = tempfile.mkstemp()
    fh = os.fdopen(fh, "w")
    SeqIO.write(seqrecords, fh, "fasta")
    fh.close()
    return fn

class XBASESeqRecord(SeqRecord):
    # a nice class to add useful functions

    def dump_to_temporary_file(self, fmt='fasta'):
        return dump_to_temporary_file((self,), fmt)

def cast_XBASESeqRecord(seqrecord):
    seqrecord.__class__ = XBASESeqRecord
    return seqrecord

def features_to_seqrecords(record):
    proteins = [get_protein(record, feature) for feature in record.features if feature.type == 'CDS' and not pseudo(feature)]
    return [cast_XBASESeqRecord(SeqRecord(Seq(p['sequence']), id=p['identifier'], description=p['description'])) for p in proteins]

def features_to_gene_seqrecords(record, offset):
    genes = [get_gene(record, feature, offset) for feature in record.features if feature.type == 'CDS']
    return [cast_XBASESeqRecord(SeqRecord(Seq(g['sequence']), id=g['identifier'], description=g['description'])) for g in genes]

class OneBasedSystem(int):
    # an integer for sequences which is 1-based (e.g. first character of string = offset 1)
    pass

class ZeroBasedSystem(int):
    pass

class WrongSystemUsed(Exception):
    pass

def find_features_at_base(seq, base, types=None):
    if not isinstance(base, OneBasedSystem):
        raise WrongSystemUsed

    features = []
    for f in seq.features:
        if (f.location.nofuzzy_start+1) <= base and f.location.nofuzzy_end >= base:
            if not types or f.type in types:
                features.append(f)
    return features

def get_snp_mutation_effect(reference_seq, base_position, base_change, gene_start, gene_end, strand):
    if not isinstance(base_position, OneBasedSystem) or not isinstance(gene_start, OneBasedSystem) or not isinstance(gene_end, OneBasedSystem):
        raise WrongSystemUsed

    ref_substr = get_seq_1_based(reference_seq, gene_start, gene_end, 1)

    # everything is 0-based from now
    relative_position = base_position - gene_start
    cns_substr = "%s%s%s" % (ref_substr[:relative_position], base_change, ref_substr[relative_position+1:])

    if strand < 0:
        ref_substr = reverse_complement(ref_substr)
        cns_substr = reverse_complement(cns_substr)

        relative_position = len(ref_substr) - relative_position - 1
    acid_position = relative_position - (relative_position % 3)
    acid_position /= 3

    dict = {}
    dict['nuc1'] = ref_substr
    dict['nuc2'] = cns_substr
#    dict['align'] = simple_alignment(ref_substr, cns_substr)

    from Bio import Translate
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    bacteria_translator = Translate.ambiguous_dna_by_id[11]

    dict['protein1'] = bacteria_translator.translate_to_stop(Seq(ref_substr, IUPAC.ambiguous_dna)).tostring() + '*'
    dict['protein2'] = bacteria_translator.translate_to_stop(Seq(cns_substr, IUPAC.ambiguous_dna)).tostring() + '*'
#    dict['protalign'] = simple_alignment(dict['protein1'], dict['protein2'])

    change = ''
    score = None

    from Bio.SubsMat.MatrixInfo import blosum62

    if 3 * len(dict['protein1']) != len(ref_substr):
        # e.g. if sequence does not encode a full-length protein, e.g. selenocysteine
        change = 'unknown'
        reference_aa = '?'
        consensus_aa = '?'
    elif len(dict['protein1']) > len(dict['protein2']):
        change = 'termination'
        reference_aa = dict['protein1'][acid_position]
        consensus_aa = '*'
    else:
        if dict['protein1'] == dict['protein2']:
            change = 'synonymous'
        else:
            change = 'non-synonymous'
        reference_aa = dict['protein1'][acid_position]
        consensus_aa = dict['protein2'][acid_position]    
        try:
            score = blosum62[(reference_aa, consensus_aa)]
        except KeyError:
            try:
                score = blosum62[(consensus_aa, reference_aa)]
            except KeyError:
                score = None

        if reference_aa == '*' and consensus_aa != '*':
            change = 'elongation'

    if change == 'non-synonymous' and reference_aa == consensus_aa:
        print "sanity error!"
        print reference_aa
        print consensus_aa
        print dict['protein1']
        print dict['protein2']
        raise SystemExit

    return change, reference_aa, consensus_aa, score

def location_sort(a, b):
    pos_cmp = cmp(a.location.nofuzzy_start, b.location.nofuzzy_start)
    if pos_cmp:
        return pos_cmp
    rank = {'gene' : 1, 'CDS' : 2}
    rank_a = rank.get(a.type, 0)
    rank_b = rank.get(b.type, 1)
    return cmp(rank_a, rank_b)

def records_to_dict(fh, type):
    hash = {}
    for rec in SeqIO.parse(fh, type):
        if rec.id in hash:
            raise ValueError('non-unique record id found')
        hash[rec.id] = rec
    return hash

