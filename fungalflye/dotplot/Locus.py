#!/usr/bin/env python3
"""This module defines Locus, the core genome coordinate class that
HistoBase is built on.
"""

from MsvUtil import prototypeable

@prototypeable
class Locus:
    """A Locus is a specific location in a genome.  It is very
    similar to a GFF feature but with fewer explicit fields.
    Public member data:
    ref = Reference sequence name (e.g. contig)
    start = Offset of first nucleotide relative to ref, starting at 1
    stop = Offset of last nucleotide
    strand = + or -
    genome = Genome object containing reference sequence
    """
    def __init__(self, ref, start, stop, strand, genome):
        """Simple constructor.  Enforces stop >= start by flipping
        strand if necessary."""
        self.ref = ref        
        self.start = int(start) 
        self.stop = int(stop)   
        self.strand = strand  
        self.genome = genome  

        if(self.start > self.stop):
            if(self.strand == "+"):
                self.strand = "-"
            else:
                self.strand = "+"
            (self.start, self.stop) = (self.stop, self.start)

    @property
    def five_prime(self):
        """Return strand-oriented first position"""
        if(self.strand == "+"):
            return start
        else:
            return stop

    @property
    def three_prime(self):
        """Return strand-oriented last position"""
        if(self.strand == "+"):
            return stop
        else:
            return start

    @classmethod
    def fromCoords(cls, coords, genome = None):
        return ParseCoords(coords, genome)

    def Locus(self):
        return self

    def Coords(self):
        """Return gbrowse formatted coordinates for this locus."""
        if(self.strand == "+"):
            return "%s:%d..%d" % (self.ref, self.start, self.stop)
        return "%s:%d..%d" % (self.ref, self.stop, self.start)

    def Sequence(self):
        """Return genomic sequence corresponding to this locus."""
        #This is somewhat redundant to Genome.GetSequence, but
        #avoids messy constructions like locus.Genome().GetSequence(locus).
        #The simple construction may make it easier to apply generic
        #algorithms to heterogenous collections of things that have sequences.
        
        return self.genome.GetSequence(self)

    # Note: This is saved for historical reference only.  As far as I
    #       know, Agilent never handled the ref field of this format
    #       correctly.
    def OldAgilentFormat(self):
        """eArray format representation: ref:5'-3'"""
        if(self.strand == "-"):
            return "%s:%d-%d" % (self.ref,
                                 self.stop,
                                 self.start)
        else: # "+" or undefined
            return "%s:%d-%d" % (self.ref,
                                 self.start,
                                 self.stop)

    def AgilentFormat(self):
        """eArray format representation: ref:5'-3'"""
        import re
        def fixref(x):
            x = re.sub(r"_",r"u",x)
            x = re.sub(r"\.",r"p",x)
            if(re.search("[\W]", x) != None):
                raise ValueError("Invalid Agilent ref: %s" % x)
            return x
            
        if(self.strand == "-"):
            return "chr%s:%d-%d" % (fixref(self.ref),
                                 self.stop,
                                 self.start)
        else: # "+" or undefined
            return "chr%s:%d-%d" % (fixref(self.ref),
                                 self.start,
                                 self.stop)

    def __repr__(self):
        # Gbrowse style representation
        return "%s:%d..%d (%s)" % (self.ref,
                                self.start,
                                self.stop,
                                self.strand)
    def __str__(self):
        return self.__repr__()

    def __len__(self):
        """Return length of locus in nucleotides."""
        return self.stop - self.start + 1

    def midpoint(self):
        """Return the center coordinate of this locus, rounded down to the
        nearest integer.
        """
        return (self.start+self.stop)//2
        
    def Overlap(self, locus):
        """Return overlap between two loci, in base pairs.

        Overlap is defined as 0 for loci on different contigs
        and <= 0 for non-overlapping loci on the same contig.
        Strand and genome are not considered."""

        if(self.ref != locus.ref):
            return 0
        
        return min(self.stop, locus.stop) - max(self.start, locus.start) + 1

    def pad(self, pad, side = None):
        """Return a copy of this Locus, padded by pad bases on both sides,
        clipped to contig boundaries (if available).
        If side is set to 5' or 3', pad only on the given side.
        """
        if((side is None) or
           (side == "5'" and self.strand == "-") or
           (side == "3'" and self.strand == "+")):
            s = self.stop+pad
            try:
                stop = min(s, len(self.genome.GetContig(self.ref)))
            except:
                stop = s
        else:
            stop = self.stop

        if((side is None) or
           (side == "5'" and self.strand == "+") or
           (side == "3'" and self.strand == "-")):
            start = max(1,self.start-pad)
        else:
            start = self.start

        return Locus.fromPrototype(
            self,
            start = start,
            stop = stop)

    def side(self, pad, side):
        """Return a locus as for Locus.pad(pad, side) but omitting the
        locus part.  (e.g., Locus.pad(500,"5'") gives a 500 bp promoter)
        """
        assert(side is not None)
        if((side == "5'" and self.strand == "-") or
           (side == "3'" and self.strand == "+")):
            s = self.stop+pad
            try:
                stop = min(s, len(self.genome.GetContig(self.ref)))
            except:
                stop = s
        else:
            stop = max(1,self.start-1)

        if((side == "5'" and self.strand == "+") or
           (side == "3'" and self.strand == "-")):
            start = max(1,self.start-pad)
        else:
            s = self.stop+1
            try:
                start = min(s, len(self.genome.GetContig(self.ref)))
            except:
                start = s

        return Locus.fromPrototype(
            self,
            start = start,
            stop = stop)
    
    def trim(self, trim, side = None):
        """Return a copy of this Locus, trimmed by trim bases on both sides.
        Raise ValueError if Locus would be entirely trimmed.
        If side is set to 5' or 3', trim only on the given side.
        """
        if(side is None):
            if(trim >= 2*len(self)):
                return None
        elif(trim >= len(self)):
            return None
        
        if((side is None) or
           (side == "5'" and self.strand == "-") or
           (side == "3'" and self.strand == "+")):
            stop = self.stop - trim
        else:
            stop = self.stop

        if((side is None) or
           (side == "5'" and self.strand == "+") or
           (side == "3'" and self.strand == "-")):
            start = self.start + trim
        else:
            start = self.start

        return Locus.fromPrototype(
            self,
            start = start,
            stop = stop)

    # TODO: decide whether to fold genome into this,
    #       add __hash__
    #       figure out if any of these can be implicitly dervied by python
    def __le__(self, rhs):
        return ((self.ref, self.start, self.stop, self.strand) <=
                (rhs.ref, rhs.start, rhs.stop, rhs.strand))
        
    def __lt__(self, rhs):
        return ((self.ref, self.start, self.stop, self.strand) <
                (rhs.ref, rhs.start, rhs.stop, rhs.strand))

    def __ge__(self, rhs):
        return ((self.ref, self.start, self.stop, self.strand) >=
                (rhs.ref, rhs.start, rhs.stop, rhs.strand))
        
    def __gt__(self, rhs):
        return ((self.ref, self.start, self.stop, self.strand) >
                (rhs.ref, rhs.start, rhs.stop, rhs.strand))

    def __ne__(self, rhs):
        return ((self.ref, self.start, self.stop, self.strand) !=
                (rhs.ref, rhs.start, rhs.stop, rhs.strand))
        
    # It would be nice to have arithmetic defined for loci so
    #   that we could make transformation to different frames
    #   of reference more obvious.  Unfortunately, there are
    #   some annoying ambiguities: e.g, if locus + int moves
    #   locus.start and locus.stop by int, should it be the
    #   same or opposite direction for + vs. - strand genes?

    #def __add__(self, rhs):
    #    """If rhs is a type where (int + rhs -> int) is defined,
    #    return a Locus with start and stop increased by rhs.
    #    (Should throw a type exception otherwise)"""
    #    return Locus(ref = self.ref,
    #                 start = self.start + rhs,
    #                 stop = self.stop + rhs,
    #                 strand = self.strand,
    #                 genome = self.genome)
    # 
    #def __sub__(self, rhs):
    #    return self.__add__(-rhs)

def unfixref(x):
    import re
    x = re.sub(r"p", r".", x)
    x = re.sub(r"u", r"_", x)
    return x

def ParseAgilent(s, genome = None):
        """Return the Locus corresponding to a string generated
        by Locus.AgilentFormat."""
        
        import re

        p = re.search("^chr(?P<ref>[^:]+):(?P<start>[\d]+)-(?P<stop>[\d]+)", s)

        # Note that Locus.__init__ should do the correct strand flipping
        # for stop > start
        return Locus(ref = unfixref(p.group("ref")),
                     start = int(p.group("start")),
                     stop = int(p.group("stop")),
                     strand = "+",
                     genome = genome)

def ParseCoords(s, genome = None):
        """Return the Locus corresponding to a string generated
        by Locus.Coords."""
        
        import re

        # Temporary hack allowing negative values
        p = re.search("^(?P<ref>[^:]+):(?P<start>-?[\d]+)\.\.(?P<stop>-?[\d]+)", s)

        # Note that Locus.__init__ should do the correct strand flipping
        # for stop > start
        return Locus(ref = p.group("ref"),
                     start = int(p.group("start")),
                     stop = int(p.group("stop")),
                     strand = "+",
                     genome = genome)

def spanningLocus(loci, genome = None, strand = "+"):
    """Return a locus spanning an iterable of loci.  Use user supplied
    strand if given, or"+" by default.  If strand is "auto", assert that all
    loci are on the strand and use that strand for the returned
    locus."""
        
    gen = iter(loci)
    locus = next(gen).Locus()
    
    if(genome is None):
        genome = locus.genome

    if(strand == "auto"):
        strand = locus.strand
        strand_flag = True
    else:
        strand_flag = False
        
    retval = Locus(ref = locus.ref,
                   start = locus.start,
                   stop = locus.stop,
                   strand = strand,
                   genome = genome)

    for i in gen:
        locus = i.Locus()
        assert(locus.ref == retval.ref)
        if(strand_flag):
            assert(locus.strand == strand)
        if(locus.start < retval.start):
            retval.start = locus.start
        if(locus.stop > retval.stop):
            retval.stop = locus.stop

    return retval

if(__name__ == "__main__"):
    locus = Locus("contig",10,1000,"+",None)
    print(locus)
