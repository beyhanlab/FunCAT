#!/usr/bin/env python3

from .Genome import Genome, Contig, MemGeneSet
from FastaFile import FastaFile
from .Sequence import DnaSequence

class FastaGenome(Genome):
    """Genome assembly represented by FASTA file.
    """
    def __init__(self, fp, name, gff3 = None):
        self.assembly = dict((name, DnaSequence(seq))
                             for (name, seq) in FastaFile(fp).seqs.items())
        self.contigs = dict((key,
                             Contig(key, len(val), assembly = self.assembly,
                                    genome = self))
                            for (key,val) in self.assembly.items())
        self.name = name

        if(gff3 is not None):
            if(isinstance(gff3,str)):
                gff3 = open(gff3)
            
            self.genes = MemGeneSet.fromGbrowse2Gff3(
                gff3, genome = self)

        self.gene_index = None
        
    def Contigs(self):
        return self.contigs.items()

if(__name__ == "__main__"):
    print("Hello, world")
