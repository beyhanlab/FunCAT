"""
Dotplot utilities used by FunCAT.

Provides genome parsing, FASTA handling, MUMmer alignment parsing,
and helper utilities for generating genome dotplots.
"""

from .FastaGenome import FastaGenome
from .FastaFile import FastaFile
from .Genome import Genome
from .GenomeCoord import GenomeCoord
from .Locus import Locus
from .Sequence import DnaSequence
from .MUMmerTools import NucmerMap, ClickCoords

__all__ = [
    "FastaGenome",
    "FastaFile",
    "Genome",
    "GenomeCoord",
    "Locus",
    "DnaSequence",
    "NucmerMap",
    "ClickCoords"
]