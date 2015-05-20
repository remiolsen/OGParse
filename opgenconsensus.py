import argparse
from intervaltree import Interval, IntervalTree
from itertools import groupby
from opgenxml import OpGenMapPlacement, OpGenRestrictionMap
from Bio import SeqIO

reverse_complements = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "Y": "R",
    "R": "Y",
    "W": "W",
    "S": "S",
    "K": "M",
    "M": "K",
    "D": "H",
    "H": "D",
    "V": "B",
    "B": "V",
    "N": "N"
}

def revcom(s):
    return complement(s)[::-1]


def complement(s):
    letters = list(s)
    complemented = []
    for base in letters:
        if base.title() in reverse_complements.keys():
            complemented.append(reverse_complements[base.title()])
        else:
            complemented.append("N")
    return ''.join(complemented)


def print_contigs(Maps, output):
    with open("{}.fasta".format(output), "w") as final_assembly:
        for Map, sequence in Maps.iteritems():
            final_assembly.write(">NouGAT_{}\n".format(Map))
            final_assembly.write("{}\n".format("".join(sequence)))


parser = argparse.ArgumentParser(
        description='generate consensus fasta files from opgen map placements')
parser.add_argument('--xml', type=str, help='Opgen XML placement file')
parser.add_argument('--fasta', type=str, help='Sequence assembly fasta file ')
parser.add_argument('--output', type=str, help='Output fasta file')
args = parser.parse_args()

opmaps = OpGenRestrictionMap(args.xml)
placements = OpGenMapPlacement(args.xml, opmaps)

opseq = {}
for omap in opmaps.notinsilico:
    mlen = opmaps.fragment_coords(omap, 0, len(opmaps.fragments[omap]))[1]
    opseq[omap] = ["n"] * mlen


fa_dict = {}
fa_file = open(args.fasta, "rU")
for record in SeqIO.parse(fa_file, "fasta"):
    fa_dict[record.id] = record.seq

to_out = []
for map_range in sorted(placements.map_fasta, reverse=True):
    m_start = map_range.begin
    m_end = map_range.end
    m_name = map_range.data.data[2]
    m_orient = map_range.data.data[3]
    f_start = map_range.data.begin
    f_end = map_range.data.end
    f_name = map_range.data.data[0]
    f_orient = map_range.data.data[1]
    f_seq = fa_dict[f_name][f_start:f_end]
    if f_orient != m_orient:
        f_seq = revcom(f_seq)
    opseq[m_name][m_start:m_end] = list(f_seq)
    to_out.append('\t'.join([m_name, str(m_start), str(m_end), f_name, str(f_start), str(f_end), f_orient]))

print "map_name\tmap_start\tmap_end\tfasta_name\tfasta_start\tfasta_end\torientation"
print "\n".join(sorted(to_out, key=lambda x: x.split("\t")[0]))
print_contigs(opseq, args.output)

    
