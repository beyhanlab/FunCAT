#!/usr/bin/env python3

# Note: LocusTransforms do not normally apply contig truncation.
#       we could add this as a post-processing step in the base
#       class, but it is more efficient to simply let it be a
#       specific subclass (this makes it easy to maintain temprorarily
#       out-of-bounds coordinates and to reserve clipping to the last
#       step in the chain).

# Note: currently using Locus.fromPrototype throughout to avoid
#       side effects on the input.  This _does_ imply a lot of
#       constructor calls for long chains.  We could try to avoid
#       this by doing a single instantiation at the begining of the
#       chain and doing in-place modifications everywhere else
#       (except for splits).  Cleanest way to make room for this
#       would be with an "in-place" flag (which could be in __call__,
#       __init__, or in a factory function mapping to two versions of
#       each class or that applies an optional "copy" decorator...

from bisect import bisect_left

from .Locus import Locus, spanningLocus

class LocusTransform:
    """A transform between two sequence coordinate spaces.
    Abstract base class."""
    def __call__(self, loci):
        """Yield iterable of transformed loci.
        Input is an iterable of loci.
        There is not necessarily a one-to-one mapping
          between input and output loci (so, if you
          need to track input/output correspondence,
          make separate calls for each input).
        """
        raise NotImplementedError

    def simple_transform(self, locus):
        """Return result of one-to-one locus transformation, or
        raise an AssertionError if the mapping is not one-to-one."""
        retval = list(self([locus]))
        assert(len(retval) == 1)
        return retval[0]

class ChainTransform(LocusTransform):
    def __init__(self, transforms = None):
        """Init with optional iterable of trnasforms."""
        if(transforms is not None):
            self.transforms = list(transforms)
        else:
            self.transforms = []

    def append(self, transform):
        """Add transform to the end of the chain."""
        self.transforms.append(transform)

    def prepend(self, transform):
        """Add transform to the beginning of the chain."""
        self.transforms = [transform]+self.transforms

    def __call__(self, loci):
        """Yield net of transform chain applied left to right."""
        for i in loci:
            L = (i,)
            for j in self.transforms:
                L = j(L)
            for j in L:
                yield j

# Note: Given that ChainTransform works correctly for the empty case,
#       there might be no need for IdentityTransform...
class IdentityTransform(LocusTransform):
    def __init__(self):
        pass
    def __call__(self, loci):
        """Yield input unmodified (for use as placeholder)."""
        for i in loci:
            yield i

class OrientTransform(LocusTransform):
    def __init__(self, reflocus, ref = None):
        self.reflocus = reflocus.Locus()
        self.ref = ref
    def __call__(self, loci):
        """Yield input flipped and shifted such that reflocus
        has coordinates 1..n in 5'->3' "+" orientation.
        (I.e., yield coordinates suitable for a GFF feature
        with reflocus as reference sequence).
        """
        for i in loci:
            L = i.Locus()
            if(self.ref is not None):
                ref = self.ref
            else:
                ref = L.ref
            if(self.reflocus.strand == "+"):
                yield Locus.fromPrototype(
                    i,
                    ref = ref,
                    start = L.start - self.reflocus.start + 1,
                    stop = L.stop - self.reflocus.start + 1)
            else: #self.reflocus.strand == "-"
                if(L.strand == "+"):
                    s = "-"
                else:
                    s = "+"
                yield Locus.fromPrototype(
                    i,
                    ref = ref,
                    # stop is the new start
                    start = self.reflocus.stop - L.stop + 1,
                    stop = self.reflocus.stop - L.start + 1,
                    strand = s)

class InverseOrientTransform(LocusTransform):
    def __init__(self, reflocus):
        self.reflocus = reflocus.Locus()
    def __call__(self, loci):
        """Yield input flipped and shifted such that 1..(len(reflocus))
        has coordinates reflocus. I.e.
        InverseOrientTransform(L)(OrientTransform(L)(x)) == x
        """
        for i in loci:
            L = i.Locus()
            if(self.reflocus.strand == "+"):
                yield Locus.fromPrototype(
                    self.reflocus,
                    start = L.start + self.reflocus.start - 1,
                    stop = L.stop + self.reflocus.start - 1,
                    strand = L.strand)
            # Note that OrientTransform with minus-strand reflocus is
            #   its own inverse
            else: #self.reflocus.strand == "-"
                if(L.strand == "+"):
                    s = "-"
                else:
                    s = "+"
                yield Locus.fromPrototype(
                    self.reflocus,
                    # stop is the new start
                    start = self.reflocus.stop - L.stop + 1,
                    stop = self.reflocus.stop - L.start + 1,
                    strand = s)
                        
class StrandedShiftTransform(LocusTransform):
    def __init__(self, shift):
        self.shift = shift
    def __call__(self, loci):
        """Yield input shifted by a 5' offset (so, + strand inputs
        shift right and - strand inputs shift left).
        """
        for i in loci:
            L = i.Locus()
            if(L.strand == "-"):
                yield Locus.fromPrototype(
                    i,
                    start = L.start - self.shift,
                    stop = L.stop - self.shift)
            else:
                # L.strand == +
                yield Locus.fromPrototype(
                    i,
                    start = L.start + self.shift,
                    stop = L.stop + self.shift)

class OrderedSegmentTransform(LocusTransform):
    """An OrderedSegmentTransform is defined by
    lists of source and target loci
    where:
    * The lists are the same length
    * The lists have a consistent, non-overlapping ordering
    * Corresponding loci have the same size

    The corresponding transform takes each segment of a query locus
    that overlaps a source locus and takes it into the corresponding
    target locus.
    """
    def __init__(self, source_loci, target_loci):
        assert(len(source_loci) == len(target_loci))
        self.source_locus = spanningLocus(source_loci)
        self.transforms = []

        #self.base2transform = [None]*len(self.source_locus)
        for (source,target) in zip(source_loci, target_loci):
            assert(len(source) == len(target))
            T = ChainTransform((OrientTransform(source),
                                InverseOrientTransform(target)))
            self.transforms.append(T)

    def __call__(self, loci):
        for locus in loci:
            for tlocus in self.transform_locus(locus):
                yield tlocus

    def dump(self):
        """Print diagnostic dump of this transform."""
        # (Not implementing this as __repr__ because it is potentially very large)
        print(self.source_locus)
        for i in self.transforms:
            print(i.transforms[0].reflocus,len(i.transforms[0].reflocus),
                  i.transforms[1].reflocus,len(i.transforms[1].reflocus))
                
    def transform_locus(self, locus):
        # Completely out of bounds loci give null transform
        if((locus.ref != self.source_locus.ref) or
           (locus.stop < self.source_locus.start) or
           (locus.start > self.source_locus.stop)):
            return []

        loci = []

        for T in self.transforms:
            source = T.transforms[0].reflocus
            
            # Note that we could speed up arriving at the first good
            # source via bisection
            if(locus.start > source.stop):
                continue
            if(locus.stop < source.start):
                break
            
            loci.append(T.simple_transform(
                Locus.fromPrototype(locus,
                                    start = max(locus.start, source.start),
                                    stop = min(locus.stop, source.stop))))

        # Note that this result is probably shuffled relative to what we
        # want for the case where locus is on the opposite strand from
        # source_locus
        return loci

class UnspliceTransform(LocusTransform):
    def __init__(self, loci):
        self.stops = []
        self.transforms = []
        s = 0
        for locus in loci:
            s += len(locus)
            self.stops.append(s)
            self.transforms.append(InverseOrientTransform(locus))
            
    def __call__(self, loci):
        for locus in loci:
            for tlocus in self.transform_locus(locus):
                yield tlocus
                
    def transform_locus(self, locus):
        # TODO: add 1-to-1 init flag to the "unsplice" family of transforms.
        #       If set, the transform would throw an exception on out-of-bounds
        #       or clipping, such that the client can trust that any successful
        #       call yields a 1-to-1 mapping at the residue level.
        #       (Alternatively, make a good interface for handling 1-to-1
        #        operations at the level of the transform classes)
        
        # Completely out of bounds loci give null transform
        if((locus.stop < 1) or (locus.start > self.stops[-1])):
            return []
        
        # Clip to transformable range
        locus = Locus.fromPrototype(locus,
                                    start = max(locus.start, 1),
                                    stop = min(locus.stop, self.stops[-1]))

        x = locus.start
        i = bisect_left(self.stops,x)
        loci = []
        while(locus.stop > self.stops[i]):
            if(i == 0):
                left = x
                right = self.stops[i]
            else:
                left = x-self.stops[i-1]
                right = self.stops[i]-self.stops[i-1]
            loci.append(self.transforms[i].simple_transform(Locus.fromPrototype(locus, start = left, stop = right)))
            x = self.stops[i] + 1
            i += 1
        if(i == 0):
            left = x
            right = locus.stop
        else:
            left = x-self.stops[i-1]
            right = locus.stop-self.stops[i-1]
        loci.append(self.transforms[i].simple_transform(Locus.fromPrototype(locus, start = left, stop = right)))
        return loci

class UnspliceCdsTransform(UnspliceTransform):
    def __init__(self, gene):
        UnspliceTransform.__init__(self, gene.CdsLoci())

class UnspliceExonsTransform(UnspliceTransform):
    def __init__(self, gene):
        UnspliceTransform.__init__(self, gene.ExonLoci())

class UnspliceProteinTransform(LocusTransform):
    def __init__(self, gene):
        self.unsplice = UnspliceCdsTransform(gene)
    def __call__(self, loci):
        return self.unsplice(
            Locus(None, (locus.start-1)*3+1, locus.stop*3, "+", None)
            for locus in loci)

class SpliceTransform(LocusTransform):
    def __init__(self, loci):
        self.locus = spanningLocus(loci, strand = "auto")
        # TODO: should do this as a deep copy, casting to Locus
        self.loci = loci[:]

        self.transforms = []
        offset = 0
        for locus in loci:
            self.transforms.append(ChainTransform(
                (OrientTransform(locus),
                 StrandedShiftTransform(offset))))
            offset += len(locus)

    @classmethod
    def invertUnspliceTransform(cls, T):
        return cls([i.reflocus for i in T.transforms])
            
    def __call__(self, loci):
        for locus in loci:
            for tlocus in self.transform_locus(locus):
                yield tlocus
                
    def transform_locus(self, locus):
        if(locus.strand != self.locus.strand):
            raise NotImplementedError

        loci = []

        for (unspliced,T) in zip(self.loci,self.transforms):
            if(unspliced.Overlap(locus) > 0):
                loci.append(T.simple_transform(Locus.fromPrototype(
                    locus,
                    start = max(locus.start, unspliced.start),
                    stop = min(locus.stop, unspliced.stop))))

        return loci

if(__name__ == "__main__"):
    print("Hello, world")
