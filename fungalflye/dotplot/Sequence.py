#!/usr/bin/env python3
"""Protein and nucleotide sequence classes optimized for string-like behaviour
and transformation (complementation, translation, etc), specialized for the
quirks of upstream fungal genome annotations (particularly CGD's C. albicans
annotation).
"""

import re
alpha_aa = "ACDEFGHIKLMNPQRSTVWY"

class Sequence:
    """A case insensitive protein or nucleotide sequence."""
    def __init__(self, seq):
        self.seq = str(seq).upper()

    def __str__(self):
        return self.seq

    def __repr__(self):
        return self.__str__()

    def __add__(self, rhs):
        """Concatenate this sequence with any object that
        can be converted to a string."""
        return Sequence(self.seq + str(rhs))

    # NOTE: At the cost of some speed, should add uppercasing
    #       on all of the casts.

    # NOTE: Comparisons are currently very permissive.  We may
    #       want type specificity
    #        e.g. DnaSequence("ATG") != ProteinSequence("ATG")
    #       On the other hand, throwing protein and dna sequences
    #       in a common grab bag is already asking for trouble.

    def __eq__(self, rhs):
        return self.seq.__eq__(str(rhs))

    def __le__(self, rhs):
        return self.seq.__le__(str(rhs))

    def __lt__(self, rhs):
        return self.seq.__lt__(str(rhs))

    def __gt__(self, rhs):
        return self.seq.__gt__(str(rhs))

    def __ge__(self, rhs):
        return self.seq.__ge__(str(rhs))

    def __ne__(self, rhs):
        return self.seq.__ne__(str(rhs))

    def count(self, sub):
        return self.seq.count(str(sub))

    def __getitem__(self, i):
        if(isinstance(i, int)):
            return self.seq[i]
        else:
            return self.__class__(self.seq[i.start:i.stop])

    def __len__(self):
        return len(self.seq)

    def find(self, sub):
        return self.seq.find(str(sub))

    def FormatFasta(self, name = "Sequence", annotation = None, w = 80):
        line1 = ">"+name
        if(annotation != None):
            line1 += " "+annotation
        line1 += "\n"
        return (line1+
                re.sub(r"(.{%d})" % w, r"\1\n", self.seq)+
                "\n")

# DNA global dicts

_universal_code = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",

    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",

    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",

    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

_ambiguity_codes = dict(
    #(i,re.compile(j)) for (i,j) in
    (("R","[AG]"),
     ("Y","[CT]"),
     ("W","[AT]"),
     ("S","[GC]"),
     ("K","[GT]"),
     ("M","[AC]"),
     ("B","[CGT]"),
     ("D","[AGT]"),
     ("H","[ACT]"),
     ("V","[ACG]"),
     ("N","[ACGT]"),
     # for completeness/convenience:
     ("A","A"),
     ("T","T"),
     ("C","C"),
     ("G","G")))

class TranslationTable:
    """{codon => amino acid} translation table with
    lazy "learning" of ambiguous codons."""
    def __init__(self, base_code):
        # make a copy of the base_code dictionary for quick lookup
        # and caching
        self.code = base_code.copy()
        # freeze a copy of the unambiguous code for resolving
        # ambiguous codons later
        self.base_codons = tuple(self.code.items())
    def __getitem__(self, codon):
        try:
            return self.code[codon]
        except KeyError:
            if(len(codon) != 3):
                return "X"
            p = re.compile("^"+"".join(_ambiguity_codes[i] for i in codon)+"$")
            aa = set(val for (key,val) in self.base_codons
                     if(p.search(key) is not None))
            if(len(aa) > 1):
                return "X"
            assert(len(aa) == 1)
            a = list(aa)[0]
            self.code[codon] = a
            return a

# TODO: Check that this doesn't break agreement with, e.g., GSC
_genetic_code = TranslationTable(_universal_code)

# C. albicans genetic code (CTG->S).  
# https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes#SG12
# Translating initiator CTG->M for SG12 is implemented in Translate below
_SG12_base = _universal_code.copy()
_SG12_base["CTG"] = "S"

_SG12 = TranslationTable(_SG12_base)

# P. tannophilus genetic code (CTG->A)
#https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes#SG26
# Allows CTG start, as for Ca
_SG26_base = _universal_code.copy()
_SG26_base["CTG"] = "A"

_SG26 = TranslationTable(_SG26_base)

# NOTE: we _could_ automatically infer these as well...
_comp = {"A":"T","T":"A","C":"G","G":"C","N":"N",
        # ambiguity codes (c.f. page 21 of the BLAST book)
        "R":"Y","Y":"R","W":"W","S":"S","K":"M","M":"K",
        "B":"V","V":"B","D":"H","H":"D"}

class DnaSequence(Sequence):
    """A single stranded DNA sequence."""

    def __init__(self, seq):
        Sequence.__init__(self, seq)

    def Complement(self):
        """Return complementary sequence.  Throws a KeyError if the
        sequence contains any non-standard letters."""
        return DnaSequence("".join(_comp[x] for x in reversed(self.seq)))

    def Antisense(self):
        """Return antisense (3'->5') sequence.  Throws a KeyError if the
        sequence contains any non-standard letters."""
        return DnaSequence("".join(_comp[x] for x in self.seq))

    def Translate(self, genetic_code = None, from_start = False):
        """Return the protein sequence obtained by translating this
        sequence with the universal genetic codes.  Ambiguous codons
        are translated as X, and stop codons as *.  The entire sequence
        is translated, so this is not as naunced as biological
        translation.  Trailing partial codons are treated as ambiguous
        codons."""

        if(genetic_code is None):
            code = _genetic_code
        elif(genetic_code == "SG12"):
            code = _SG12
        else:
            raise ValueError("Unknown genetic code %s" % genetic_code)

        # Special case for non-canonical initiators.
        # We assume that the calling function knows the translation
        # initiation site.
        (first, init) = (0, "")
        if(from_start and genetic_code == "SG12"
           and self.seq.startswith("CTG")):
            (first, init) = (3, "M")

        return ProteinSequence(
            init+"".join(code[self.seq[i:i+3]]
                         for i in range(first, len(self.seq), 3)))

    def Format3frame(self, w = 80, genetic_code = None):
        """Return 3 frame translation formatted a la DNA Strider."""
        nt = str(self)
        f0 = "".join(map(lambda x: x+"  ",
                         str(self.Translate(genetic_code = genetic_code))))
        f1 = "".join(map(lambda x: " "+x+" ",
                         str(self[1:].Translate(genetic_code = genetic_code))))
        f2 = "".join(map(lambda x: "  "+x,
                         str(self[2:].Translate(genetic_code = genetic_code))))
        retval = ""
        i = 0
        while(i < len(nt)):
            retval += nt[i:i+w]+"\n"
            retval += f0[i:i+w]+"\n"
            retval += f1[i:i+w]+"\n"
            retval += f2[i:i+w]+"\n"
            i += w
        return retval

    # Tm member functions.  We may want to factor these into an
    # energy class EGAD/AOS style

    def GC(self):
        """Returns fraction GC = (Ng+Nc)/N"""
        return (self.seq.count("G")+self.seq.count("C"))/len(self.seq)

    def Tm1(self):
        """Calculate Tm from "the standard formula" as given in
        GenomeRes_16_271.  Assumes pure [ATGC] composition."""
        
        return 64.9 + 41.0*(self.seq.count("G")+self.seq.count("C")-16.4)/len(self.seq)

    def disjoint_cds(self, minlen = 3, startcodons = None, closecds = False,
                     bothstrands = False):
        """Given a DnaSequence object, return a list of (start,stop)
        coordinates for all open reading frames that do not overlap in a
        common frame of reference.  (start,stop) follow the convention of
        Gene.Cds(): 1-based coordinates relative to the transcript, stop
        is the last base, stop codon not included.  If startcodons is None
        then all non-stop codons are permitted.  If closecds=True,
        then final, unterminated CDS are allowed.  If bothstrands=True,
        then all six reading frames are considered (with antisense CDS
        having start>stop)."""
        cds = []

        gen = [(str(self), lambda x: x)]
        if(bothstrands):
            gen.append((str(self.Complement()), lambda x: len(self)-x+1))

        for (s,t) in gen:
            starts = [None,None,None]
            stopcodons = ("TAA","TGA","TAG")
            for i in range(len(s)-2):
                f = i % 3
                c = s[i:i+3]
                if(starts[f] is not None):
                    if(c in stopcodons):
                        # string indexing
                        if(i+3-starts[f] >= minlen):
                            cds.append((t(starts[f]+1), t(i)))
                        starts[f] = None
                elif((not c in stopcodons) and (
                    (startcodons is None) or (c in startcodons))):
                    starts[f] = i
            if(closecds):
                for j in starts:
                    if(j is not None):
                        # Note that this allows the final codon to be ragged
                        cds.append((t(j+1), t(len(s))))
        return cds

class ProteinSequence(Sequence):
    """An amino acid sequence."""
    # Molecular weights in daltons (from Creighton, page 4)
    __mass = {
        "A":71.09,
        "R":156.19,
        "N":114.11,
        "D":115.09,
        "C":103.15,
        "Q":128.14,
        "E":129.12,
        "G":57.05,
        "H":137.14,
        "I":113.16,
        "L":113.16,
        "K":128.17,
        "M":131.19,
        "F":147.18,
        "P":97.12,
        "S":87.08,
        "T":101.11,
        "W":186.21,
        "Y":163.18,
        "V":99.14,
        # Use average weight for unknown amino acid
        "X":119.40,
        # Stop codons weigh nothing =)
        "*":0.0}
    def __init__(self, seq):
        Sequence.__init__(self, seq)

    def mass(self):
        """Return molecular weight in daltons
        (not a perfect match to ExPASy's ProtParam)"""
        return sum(ProteinSequence.__mass[i] for i in self)

# Unit test
if(__name__ == "__main__"):
    dna = DnaSequence("CATATGGAGGAGTGACAT")
    import sys
    sys.stdout.write(dna.FormatFasta(name="dna",annotation="test"))
    sys.stdout.write(dna.Complement().FormatFasta(
        name="dna",annotation="complemented"))
    sys.stdout.write(dna.Antisense().FormatFasta(
        name="dna",annotation="antisense"))
    sys.stdout.write(dna.Format3frame())
    protein = dna.Translate()
    sys.stdout.write(protein.FormatFasta(name="protein",annotation="translated"))
    print("mass",protein.mass())
    print(dna.disjoint_cds())

                     
def disjoint_cds(seq, **kw):
    """Deprecated.  Use DnaSequence.disjoint_cds instead (identical behavior)"""
    return seq.disjoint_cds(**kw)


        
