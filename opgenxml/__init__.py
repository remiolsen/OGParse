import xml.etree.ElementTree as ET
from intervaltree import Interval, IntervalTree



class OpGenRestrictionMap():

    def __init__(self, xmlfile):
        self.xmlroot = ET.parse(xmlfile).getroot()
        self.fragments = {}
        self.notinsilico = []
        self._get_fragments()

    def _get_fragments(self):
        for rmap in self.xmlroot.iter('RESTRICTION_MAP'):
            mapID = rmap.get('ID').strip()
            self.fragments[mapID] = [frag.attrib for frag in rmap.iter('F')]
            if rmap.get('INSILICO') == 'false':
                self.notinsilico.append(mapID)


    def fragment_coords(self, ID, start, end):
        fragments = self.fragments[ID]
        fasta_start = sum([int(frag['S']) for frag in fragments[0:start]])
        fasta_end = sum([int(frag['S']) for frag in fragments[0:end]])

        return(min(fasta_start,fasta_end),max(fasta_start,fasta_end))


class OpGenMapPlacement():
    
    def __init__(self, xmlfile, rs_map):
        self.xmlroot = ET.parse(xmlfile).getroot()
        self.align_chunks = {}
        self.map_fasta = IntervalTree()
        self.fasta_map = IntervalTree()
        self._get_chunks(rs_map)

    def _get_chunks(self, rs_map):
        for align in self.xmlroot.findall('MAP_ALIGNMENT'):
            S1 = align.attrib['MAP1'].split()[0]
            S2 = align.attrib['MAP2'].split()[0]

            chunks = [chunk.attrib for chunk in align.findall('CHUNK')]
            start_S1 = int(chunks[0]["S1"])
            start_S2 = int(chunks[0]["S2"])
            end_S1 = int(chunks[-1]["S1"])
            end_S2 = int(chunks[-1]["S2"])

            coords_S1 = rs_map.fragment_coords(S1, start_S1, end_S1)
            coords_S2 = rs_map.fragment_coords(S2, start_S2, end_S2)
            if not self.align_chunks.has_key((S1,S2)):
                self.align_chunks[(S1,S2)] = {coords_S1: chunks}
            else:
                self.align_chunks[(S1,S2)][coords_S1] = chunks
            
            if start_S1 > end_S1:
                orientationS1 = "-"
            else:
                orientationS1 = "+"
            if start_S2 > end_S2:
                orientationS2 = "-"
            else:
                orientationS2 = "+"

            self.map_fasta[coords_S1[0]:coords_S1[1]] = Interval(coords_S2[0], coords_S2[1], (S2, orientationS2, S1, orientationS1))
            self.fasta_map[coords_S2[0]:coords_S2[1]] = Interval(coords_S1[0], coords_S1[1], (S2, orientationS2, S1, orientationS1))

