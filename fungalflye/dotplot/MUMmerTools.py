"""Utilities for parsing and viewing MUMmer output."""

from .Locus import Locus
from .GenomeCoord import GenomeCoord
from PCA import pca_plot

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import re

class NucmerMap:
    """nucmer-based map between two genomes.
    (Currently restricted to show-coords output)."""
    def __init__(self, pairs, genomeA, genomeB, coordA, coordB):
        self.genomeA = genomeA
        self.genomeB = genomeB
        self.pairs = pairs
        if(coordA is None):
            self.coordA = GenomeCoord(genomeA)
        else:
            self.coordA = coordA
        if(coordB is None):
            self.coordB = GenomeCoord(genomeB)
        else:
            self.coordB = coordB
    
    @classmethod
    def from_coords(cls, fp, genomeA, genomeB,
                    coordA = None, coordB = None):
        """Init from show-coords output."""
        # Tested on output of show-coords -r -c -l A_vs_B.nucmer.delta
        if(isinstance(fp,str)):
            fp = open(fp,"rt")
        # Burn header
        line = next(fp)
        n = 1
        while(not line.startswith("==")):
            line = next(fp)
            n += 1
        #print "Burned %d header lines" % n
        
        pairs = []
        for line in fp:
            (cA,cB,aligned_lengths,identity,total_lengths,coverage,refs) = line.split("|")
            (strandA,strandB,refA,refB) = [i.strip() for i in refs.strip().split()]
            (startA,stopA) = [int(i) for i in cA.strip().split()]
            (startB,stopB) = [int(i) for i in cB.strip().split()]
            # N.B. show-coords -d correctly flips (start,stop) for - strand
            pairs.append((Locus(refA,startA,stopA,"+",genomeA),
                          Locus(refB,startB,stopB,"+",genomeB)))
                        
        return cls(pairs,genomeA,genomeB,coordA,coordB)

    @pca_plot
    def plot(self, ax = None, transpose = False, gridA = True, gridB = True):
        coordA = self.coordA
        coordB = self.coordB
        pc = self.pair_coords
        if(transpose):
            (coordA,coordB) = (coordB,coordA)
            def pc(i):
                (x0,x1,y0,y1) = self.pair_coords(i)
                return (y0,y1,x0,x1)

        if(gridA):
            for i in coordA.offsets.values():
                ax.axvline(i, color = "purple")
        if(gridB):
            for i in coordB.offsets.values():
                ax.axhline(i, color = "purple")

        lines = []
        for i in self.pairs:
            (x0,x1,y0,y1) = pc(i)
            lines.append(((x0,y0),(x1,y1)))
        ax.add_collection(LineCollection(lines))

        ax.set_xlim(1,coordA.max)
        ax.set_ylim(1,coordB.max)

    def pair_coords(self, pair):
        (a,b) = pair
        A = self.coordA(a)
        if(A.strand == "+"):
            x0 = A.start
            x1 = A.stop
        else:
            x0 = A.stop
            x1 = A.start
        B = self.coordB(b)
        if(B.strand == "+"):
            y0 = B.start
            y1 = B.stop
        else:
            y0 = B.stop
            y1 = B.start

        return (x0,x1,y0,y1)

# Telomere related calculations and plots

from MsvUtil import median, hdict, argmax
def nucmer_stats(comp):
    sums = {}
    for (i,c) in comp.genomeA.Contigs():
        sums[i] = {}
        for (j,c) in comp.genomeB.Contigs():
            sums[i][j] = 0
    for (t,g) in comp.pairs:
        sums[t.ref][g.ref] += len(t)
    best = dict((i, argmax(sums[i].items(), lambda x: x[1])) for i in sums)
    contig2best = {}
    for (t, (g,l)) in best.items():
        try:
            contig2best[g].append(t)
        except KeyError:
            contig2best[g] = [t]
    pi = hdict(comp.pairs, lambda x: (x[0].ref,x[1].ref))
    return sums, best, contig2best, pi

def get_offset(tl,gl):
    """Given a pair of dotplot loci (with start <= stop and gl.strand="-" for antiparallel orientation)
    return offset,orientation specifying the 5' end of the full tl sequence relative to gl."""
    if(gl.strand == "+"):
        return gl.start - tl.start,"+"
    else:
        return gl.stop + tl.start,"-"

# call as estimate_offset(pi[(t,g)])
def oriented_offset(pairs):
    """Given all nucmer alignments between t and g, return the 
    median expected offset and orientation of the 5' end of t relative to g.
    """
    offsets = {"+":[],"-":[]}
    for (tl,gl) in pairs:
        offset,strand = get_offset(tl,gl)
        offsets[strand].append(offset)
        
    if(len(offsets["+"]) >= len(offsets["-"])):
        return median(offsets["+"]),"+"
    return median(offsets["-"]),"-"

def estimate_offset(pairs):
    (offset, orientation) = oriented_offset(pairs)
    return offset

def plot_tread_offsets(comp, contig2best, pi, g, merge = True):
    offsets = {"GGGTTA":[],"TAACCC":[]}
    for t in contig2best[g]:
        r = t.split("_")[-1]
        offsets[r].append(estimate_offset(pi[(t,g)]))

    (fig,ax) = plt.subplots()
    if(merge):
        ax.plot(sorted(offsets["GGGTTA"]+offsets["TAACCC"]),"bo")
    else:
        ax.plot(sorted(offsets["GGGTTA"]),"bo")
        ax.plot(sorted(offsets["TAACCC"]),"ro")
    ax.axhline(0)
    ax.axhline(len(comp.genomeB.GetContig(g)))
    return fig

import re
tel5_re = re.compile("TAACCC")
tel3_re = re.compile("GGGTTA")

class TelomerePlot:
    def __init__(self, comp, pi, g, xlim = 100000, maxgap = 100,
                 _3prime = False, contig2best = None):
        self.fig = plt.figure()
        # All events received by this frame, first to last
        self.eventbin = []
        # All events resolved to features, first to last
        self.clickbin = []
        # Map from clickable matplotlib artists to corresponding features
        self.artist2feature = {}

        y = 1
        if(_3prime):
            L = len(comp.genomeB.GetContig(g))
            print(L)
        for (t,contig) in comp.genomeA.Contigs():
            try:
                if(contig2best is not None):
                    if(t not in contig2best[g]):
                        continue
                pairs = pi[t,g]
            except KeyError:
                continue
            (o1,orientation) = oriented_offset(pairs)
            tr = comp.genomeA.GetContig(t)
            tlen = len(tr.Locus())
            if(orientation == "+"):
                offset = o1
            else:
                offset = o1 - tlen+1
            if(((not _3prime) and (offset < xlim)) or
               (_3prime and (offset >= L - xlim))):
                if(orientation == "+"):
                    fmt = "b-"
                    seq = str(tr.Locus().Sequence())
                else:
                    fmt = "r-"
                    seq = str(tr.Locus().Sequence().Complement())
                for i in plt.plot([offset,offset+len(tr)-1],[y,y],fmt,
                                  picker = True):
                    self.artist2feature[i] = tr
                for (p,c) in ((tel5_re,"c-"),(tel3_re,"m-")):
                    for i in p.finditer(seq):
                        x = i.span()[0]
                        plt.plot([offset+x,offset+x],[y-.2,y+.2],c)
                    for (tl,gl) in pairs:
                        o2,s2 = get_offset(tl,gl)
                        if((s2 == orientation) and (abs(o2-o1) <= maxgap)):
                            if(orientation == "+"):
                                plt.plot([offset+tl.start-1,offset+tl.stop-1],
                                         [y-.1,y-.1],"k-")
                            else:
                                plt.plot([offset+tlen - tl.start+1,offset+tlen-tl.stop+1],
                                         [y-.1,y-.1],"k-")

                y += 1
        y = 0
        if(_3prime):
            plt.axvline(L,color="green")
            offset = max(0,L-xlim)
            seq = str(comp.genomeB.GetContig(g)[offset:])
        else:
            plt.axvline(1,color="green")
            offset = 1
            seq = str(comp.genomeB.GetContig(g)[:xlim])
        for (p,c) in ((tel5_re,"c-"),(tel3_re,"m-")):
            for i in p.finditer(seq):
                x = i.span()[0]
                plt.plot([offset+x,offset+x],[y-.2,y+.2],c)

        self.fig.canvas.mpl_connect('pick_event', self.onpick)

    def onpick(self, event):
        self.eventbin.append(event)
        try:
            feature = self.artist2feature[event.artist]
            self.clickbin.append(feature)
        except:
            self.clickbin.append(None)

# Utility for adding interactive clicks to NucmerMap plots
# TODO: Should bundle ClickCoords with the plot in a Frame object
#       as in HistoPlot2
class ClickCoords:
    def __init__(self, comp):
        self.comp = comp
        self.eventbin = []
        self.clickbin = []
    def __call__(self, event):
        self.eventbin.append(event)
        try:
            self.clickbin.append((self.comp.coordA.invert_coord(event.xdata), 
                                  self.comp.coordB.invert_coord(event.ydata)))
        except:
            self.clickbin.append(None)

class NucmerAlignmentBlock:
    def __init__(self, locusA, locusB):
        self.locusA = locusA
        self.locusB = locusB
    def show(self, pad = 5):
        raise NotImplementedError
    def showcomp(self, pad = 5):
        raise NotImplementedError
    def seqA(self):
        raise NotImplementedError
    def seqB(self):
        raise NotImplementedError
        
class SnpAlignmentBlock(NucmerAlignmentBlock):
    def show(self, pad = 5):
        return "".join(("%s[%s]%s\n" % (
            i.side(pad, side="5'").Sequence(),
            i.Sequence(),
            i.side(pad, side="3'").Sequence())
        ) for i in (self.locusA, self.locusB))
    def showcomp(self, pad = 5):
        return "".join(("%s[%s]%s\n" % (
            i.side(pad, side="3'").Sequence().Complement(),
            i.Sequence().Complement(),
            i.side(pad, side="5'").Sequence().Complement())
        ) for i in (self.locusA, self.locusB))
    def seqA(self):
        return str(self.locusA.Sequence())
    def seqB(self):
        return str(self.locusB.Sequence())
    
class InsAlignmentBlock(NucmerAlignmentBlock):
    def show(self, pad = 5):
        return "%s[%s]%s\n%s[%s]%s\n" % (
            self.locusA.side(pad, side="5'").Sequence(),
            self.locusA.Sequence(),
            self.locusA.side(pad, side="3'").Sequence(),
            # This next bit might need to be strand-specific
            self.locusB.pad(pad-1,side="5'").Sequence(),
            " "*len(self.locusA),
            self.locusB.side(pad,side="3'").Sequence())
    def showcomp(self, pad = 5):
        return "%s[%s]%s\n%s[%s]%s\n" % (
            self.locusA.side(pad, side="3'").Sequence().Complement(),
            self.locusA.Sequence().Complement(),
            self.locusA.side(pad, side="5'").Sequence().Complement(),
            # This next bit might need to be strand-specific
            self.locusB.pad(pad-1,side="3'").Sequence().Complement(),
            " "*len(self.locusA),
            self.locusB.side(pad,side="5'").Sequence().Complement())
    def seqA(self):
        return str(self.locusA.Sequence())
    def seqB(self):
        return " "*len(self.locusA)
    
class DelAlignmentBlock(NucmerAlignmentBlock):
    def show(self, pad = 5):
        return "%s[%s]%s\n%s[%s]%s\n" % (
            self.locusA.pad(pad-1,side="5'").Sequence(),
            " "*len(self.locusB),
            self.locusA.side(pad,side="3'").Sequence(),
            self.locusB.side(pad, side="5'").Sequence(),
            self.locusB.Sequence(),
            self.locusB.side(pad, side="3'").Sequence())
    def showcomp(self, pad = 5):
        return "%s[%s]%s\n%s[%s]%s\n" % (
            self.locusA.pad(pad-1,side="3'").Sequence().Complement(),
            " "*len(self.locusB),
            self.locusA.side(pad,side="5'").Sequence().Complement(),
            self.locusB.side(pad, side="3'").Sequence().Complement(),
            self.locusB.Sequence().Complement(),
            self.locusB.side(pad, side="5'").Sequence().Complement())
    def seqA(self):
        return " "*len(self.locusB)
    def seqB(self):
        return str(self.locusB.Sequence())

class NucmerSubAlignment:
    def __init__(self, locusA, locusB, diffs):
        self.locusA = locusA
        self.locusB = locusB
        self.diffs = diffs
        self.refA = self.locusA.ref
        self.refB = self.locusB.ref
        self.genomeA = self.locusA.genome
        self.genomeB = self.locusB.genome
        
    def diff2loci(self, diff):
        # N.B.: have to track strand as it is not memorized by single base locations (which are common)
        ((startA,stopA,strandA),(startB,stopB,strandB)) = diff
        if(stopA is None):
            assert(stopB is not None)
            return DelAlignmentBlock(
                Locus(self.refA,startA,startA,strandA,self.genomeA),
                Locus(self.refB,startB,stopB,strandB,self.genomeB))
        elif(stopB is None):
            return InsAlignmentBlock(
                Locus(self.refA,startA,stopA,strandA,self.genomeA),
                Locus(self.refB,startB,startB,strandB,self.genomeB))
        return SnpAlignmentBlock(
            Locus(self.refA,startA,stopA,strandA,self.genomeA),
            Locus(self.refB,startB,stopB,strandB,self.genomeB))

    def __iter__(self):
        for i in self.diffs:
            yield self.diff2loci(i)
    
class NucmerAlignment:
    """NucmerAlignment is isomorphic to show-aligns output."""
    def __init__(self, subalignments, genomeA, refA, genomeB, refB):
        self.subalignments = subalignments
        self.genomeA = genomeA
        self.refA = refA
        self.genomeB = genomeB
        self.refB = refB

    def __iter__(self):
        for i in self.subalignments:
            for j in i:
                yield j
        
    @classmethod
    def from_show_aligns(cls, fp, genomeA, genomeB):
        retval = []
        diffs = []
        # Burn header
        line = ""
        while(not(line.startswith("="))):
            line = next(fp)

        alignment_re = re.compile(r"^-- Alignments between (?P<refA>[\S]+) and (?P<refB>[\S]+)")
        begin_re = re.compile(r"^-- BEGIN alignment \[(?P<coords>.*)\][\s]*$")
        end_re = re.compile(r"^--   END alignment \[(?P<coords>.*)\][\s]*$")
        seq_re = re.compile(r"^(?P<pos>[\d]+)(?P<space>[\s]*)(?P<seq>[\S]+)")
        var_re = re.compile(r"(?P<var>\^+)")

        def getpos(locus, consumed):
            if(locus.strand == "+"):
                return locus.start+consumed
            return locus.stop-consumed
        
        def make_diff(
                # Reference sequences
                locusA, locusB,
                # Current positions in the reference sequences
                consumedA, consumedB,
                # Start and end offsets relative to the current positions
                diff_startA, diff_stopA, diff_startB, diff_stopB,
                # Type of diff to create
                state):
            if(state in ("SNP","INS")):
                startA = getpos(locusA, consumedA + diff_startA)
                stopA = getpos(locusA, consumedA + diff_stopA)
                A = sorted((startA, stopA))+[locusA.strand]
            else:
                A = (getpos(locusA, consumedA + diff_startA - 1),
                     None,locusA.strand)

            if(state in ("SNP","DEL")):
                startB = getpos(locusB, consumedB + diff_startB)
                stopB = getpos(locusB, consumedB + diff_stopB)
                B = sorted((startB, stopB))+[locusB.strand]
            else:
                B = (getpos(locusB, consumedB + diff_startB - 1),
                     None,locusB.strand)

            return (A,B)
        
        refA = None
        refB = None
        
        for line in fp:
            p = alignment_re.search(line)
            if(p is None):
                assert(len(line.strip()) == 0)
                continue
            if(refA is None):
                refA = p.group("refA")
                refB = p.group("refB")
            else:
                assert(refA == p.group("refA"))
                assert(refB == p.group("refB"))
                
            for line in fp:
                if(line.startswith("=")):
                    break
                p = begin_re.search(line)
                if(p is None):
                    assert(len(line.strip()) == 0)
                    continue
                coords = p.group("coords")
                ((strandA,startA,hA,stopA),
                 (strandB,startB,hB,stopB)) = [i.split()
                                               for i in coords.split("|")]

                assert(hA == hB == "-")

                # N.B.: minus strand coords are already swapped, so
                # don't double-correct
                locusA = Locus(ref=refA,start=int(startA),stop=int(stopA),
                               strand="+",genome=genomeA)
                locusB = Locus(ref=refB,start=int(startB),stop=int(stopB),
                               strand="+",genome=genomeB)
                consumedA = 0
                consumedB = 0
                
                for line in fp:
                    if(len(line.strip()) == 0):
                        continue
                    p = end_re.search(line)
                    if(p is not None):
                        assert(p.group("coords") == coords)
                        retval.append(NucmerSubAlignment(locusA, locusB, diffs))
                        diffs = []
                        break
                    p = seq_re.search(line)
                    if(p is None):
                        print(line)
                        raise ValueError
                    posA = int(p.group("pos"))
                    assert(posA == getpos(locusA, consumedA))
                      
                    offset = len(p.group("pos"))+len(p.group("space"))
                    seqA = p.group("seq")
                    
                    line = next(fp)
                    p = seq_re.search(line)
                    posB = int(p.group("pos"))
                    assert(posB == getpos(locusB, consumedB))
                    assert(offset == len(p.group("pos"))+len(p.group("space")))
                    seqB = p.group("seq")
                    assert(len(seqB) == len(seqA))

                    line = next(fp)
                    j = 0
                    for var in var_re.finditer(line[offset:]):
                        old_j = j
                        (i,j) = var.span()
                        # consume the non-variant positions we just skipped
                        consumedA += i-old_j
                        consumedB += i-old_j
                        
                        # consume the current variant block
                        startA = getpos(locusA,consumedA)
                        startB = getpos(locusB,consumedB)

                        varA = seqA[i:j]
                        varB = seqB[i:j]

                        state = None
                        last_n = 0
                        deltaA = 0
                        deltaB = 0
                        dA = 0
                        dB = 0
                        for (n,(a,b)) in enumerate(zip(varA,varB)):
                            if(a == "."):
                                curstate = "DEL"
                                dA += 1
                            elif(b == "."):
                                curstate = "INS"
                                dB += 1
                            else:
                                curstate = "SNP"
                            if(curstate != state):
                                if(state is not None):
                                    diffs.append(make_diff(
                                        locusA,locusB,consumedA,consumedB,
                                        last_n+deltaA, n-1+deltaA,
                                        last_n+deltaB, n-1+deltaB,
                                        state))
                                    deltaA -= dA
                                    deltaB -= dB
                                    dA = 0
                                    dB = 0

                                state = curstate
                                last_n = n

                        if(state is not None):
                            diffs.append(make_diff(
                                locusA,locusB,consumedA,consumedB,
                                last_n+deltaA, len(varA)-1+deltaA,
                                last_n+deltaB, len(varB)-1+deltaB,
                                state))        
                                
                        consumedA += len(varA)-varA.count(".")
                        consumedB += len(varB)-varB.count(".")
                    rest = len(seqA)-j
                    assert(rest >= 0)
                    consumedA += rest
                    consumedB += rest
        
        return cls(retval, genomeA, refA, genomeB, refB)
