from Bio import SeqIO
import re
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compile a BED file of all gaps and ambiguous bases')
    parser.add_argument('fastafile', type=str, help='Fasta file')
    parser.add_argument('--ambiguous', action='store_true', default=False, 
            help='Outputs ambiguous bases. (Default false)')
    args = parser.parse_args()
    for sequence in SeqIO.parse(args.fastafile, "fasta"):
        sequence_name = sequence.id
        seqstr = str(sequence.seq)
        
        if args.ambiguous:
            refilter = re.compile('(?P<gap>[Nn]+)|(?P<amb>[YRWSKMDHVByrwskmdhvb]+)')
        else:
            refilter = re.compile('(?P<gap>[Nn]+)')

        for event in re.finditer(refilter,seqstr):
            if event.group("gap") is not None:
                etype = "gap"
            elif event.group("amb") is not None:
                etype = "ambiguous"
            print "%s\t%s\t%s\t%s" % (sequence_name, event.start(), event.end(), etype)

if __name__ == "__main__":
    main()
