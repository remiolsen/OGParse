import argparse
from intervaltree import Interval, IntervalTree
from itertools import groupby
from opgenxml import OpGenMapPlacement, OpGenRestrictionMap        


parser = argparse.ArgumentParser(
        description='Parse OpGen XML files and parse opticalmap-insilico digested fasta alignments')
parser.add_argument('--xml', type=str, help='Opgen XML placement file')
parser.add_argument('--bed_in', type=str, help='Input bed file to be translated to opgen alignment coordinates')
parser.add_argument('--bed', type=str, help='Output bed file (default: print to stdout)')
args = parser.parse_args()

opmaps = OpGenRestrictionMap(args.xml)
placements = OpGenMapPlacement(args.xml, opmaps)
bed_out = []

#Target (opgen map) specific information
tname = placements.align_chunks.keys()[0][0]
tfrags = opmaps.fragments[tname]
tchunks = (int(tfrags[0]["I"]), int(tfrags[-1]["I"]))
tsize = opmaps.fragment_coords(tname, min(tchunks), max(tchunks))[1]

#Annotate regions of fasta file as placed on optical maps
if args.bed_in is None:
    for interval in sorted(placements.map_fasta.items()):
        tstart = interval.begin
        tend = interval.end

        #Query (in-silico digested sequence) specific information
        q = interval.data
        qname = q.data[0]
        qstrand = q.data[1]
        bedline = "%s\t%s\t%s\t%s\t%s\t%s" % (tname, tstart, tend, qname, 999, qstrand)
        bed_out.append(bedline)
# We are converting bed coordinates to optical map coordintes
else:
    bedintervals = IntervalTree()
    try:
        with open(args.bed_in) as bi:
            for bi_line in bi:
                bi_line = bi_line.split()
                bedintervals[int(bi_line[1]):int(bi_line[2])] = (bi_line[0], bi_line[3])
    except IOError as e:
        print "Something wrong with input bed file: %s" % (e)
        exit()
    #bin the contigs by map coordinates
    fmap = {}
    tree = sorted(placements.fasta_map, key=lambda x: x[2])
    for k, g in groupby(tree, key=lambda x: x[2]):
        iv = list(g)
        if fmap.has_key(k.data):
            fmap[k.data].update(iv)
        else:
            fmap[k.data] = IntervalTree(iv)
    #Bin the bed coordiantes
    subbeds = {}
    btree = sorted(bedintervals, key=lambda x: x[2])
    for k, g in groupby(btree, key=lambda x: x[2]):
        if subbeds.has_key(k[0]):
            subbeds[k[0]].update(list(g))
        else:
            subbeds[k[0]] = IntervalTree(list(g))

    #Find features that overlap with placed contigs
    for fasta_key, finterval in fmap.iteritems():

        for bedline in subbeds[fasta_key[0]]:
            feat_start = bedline.begin
            feat_end = bedline.end
            overlaps = finterval[feat_start:feat_end]

            # We found overlaps, the feature need to be translated to map-coordinates
            if overlaps:
                ctg = list(overlaps)[0]
                fa_start = ctg.begin
                fa_end = ctg.end
                map_start = ctg.data.begin
                map_end = ctg.data.end
                if fasta_key[1] == "+":
                    feat_map_start = map_start + feat_start - fa_start
                    feat_map_end = map_start + feat_end - fa_start
                else:
                    feat_map_end = map_end - (feat_start - fa_start)
                    feat_map_start = map_end - (feat_end - fa_start)

                bed_line =  "%s\t%s\t%s\t%s\t%s\t%s" % (
                        tname, feat_map_start, feat_map_end, bedline.data[1], 999, fasta_key[1])
                bed_out.append(bed_line)

# We're done, writing to file or to stdout
bed_out = sorted(bed_out, key=lambda x: int(x.split()[1]))
if args.bed is None:
    for out_line in bed_out:
        print out_line
else:
    with open(args.bed, "w") as bout:
        bout.write("\n".join(bed_out))

