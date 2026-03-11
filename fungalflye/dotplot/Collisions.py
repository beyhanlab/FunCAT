#!/usr/bin/env python3
"""
Functions for 1D feature collisions
"""

from .Locus import Locus

class Collisions:
    def __init__(self, A, B):
        """Return collisions between sequence A and sequence B where
        elements x support x.start and x.stop and are sorted
        on x.start"""

        self.overlaps = []
        self.a_index = {}
        self.b_index = {}

        if((len(A) > 0) and (len(B) > 0)):
            self.__collide(A, B)

    def __overlap(self, ap, bp):
        """Add an overlap pair and update any relevant indices"""
        self.overlaps.append((ap, bp))
        try:
            self.a_index[ap].append(bp)
        except:
            self.a_index[ap] = [bp]
        try:
            self.b_index[bp].append(ap)
        except:
            self.b_index[bp] = [ap]

    def __collide(self, A, B):
        """Actual collision test"""

        bstart = ap = bp = 0

        while(ap < len(A)):
            first = None
            a = A[ap]
            b = B[bp]

            while(a.stop >= b.start):
                if(b.stop >= a.start):
                    self.__overlap(ap, bp)
                    if(first == None):
                        first = bp
                bp += 1
                if(bp >= len(B)):
                    break
                b = B[bp]

            ap += 1
            if(first != None):
                bstart = first
            bp = bstart

    def __getitem__(self, i):
        """Return the ith overlap as a tuple of (A index, B index)"""
        return self.overlaps[i]

    def __getslice__(self, i, j):
        return self.overlaps[i:j]

    def __len__(self):
        return len(self.overlaps)

    def OverlapA(self, i):
        """Return a list of indices into B for each member of B
        overlapping A[i]."""
        try:
            return self.a_index[i]
        except:
            return []

    def OverlapB(self, i):
        """Return a list of indices into A for each member of A
        overlapping B[i]."""
        try:
            return self.b_index[i]
        except:
            return []

def SortLoci(loci):
    """Given an iterable of objects supporting Locus(), return a dict of
    lists where the keys are refs and the values are lists of the input
    objects sorted on start position."""
    refs = {}
    for i in loci:
        try:
            refs[i.Locus().ref].append(i)
        except:
            refs[i.Locus().ref] = [i]
            
    for (key, val) in refs.items():
            val.sort(key = lambda x: x.Locus().start)
    return refs

class RefCollisions:
    """RefCollisions is a new implementation of Collisions.
    The changes that break (or threaten) backwards compatibility are:
    1) The initializing objects A & B can be either dicts or iterables.
       If they are dicts, the keys are assumed to be refs, different refs
       are assumed to be orthogonal, and the values are assumed to be
       lists sorted on start.
       If they are iterables, they will be sorted into appropriate dicts of
       lists.
       (A and B may be of different types from each other).
    2) Overlaps are tracked in terms of object references rather than indices.
    3) The elements are assumed to support .Locus(), but are not required
       to be of class Locus.
    """
    def __init__(self, A, B):
        """Return collisions between sequence A and sequence B where
        elements x support x.start and x.stop and are sorted
        on x.start"""

        self.__overlaps = []
        self.a_index = {}
        self.b_index = {}

        if(type(A) == dict):
            self.__a = A
        else:
            self.__a = SortLoci(A)

        if(type(B) == dict):
            self.__b = B
        else:
            self.__b = SortLoci(B)

        if((len(self.__a) > 0) and (len(self.__b) > 0)):
            for ref in frozenset(self.__a.keys()).intersection(
                frozenset(self.__b.keys())):
                
                self.__collide(self.__a[ref], self.__b[ref])

    def __overlap(self, ap, bp):
        """Add an overlap pair and update any relevant indices"""
        self.__overlaps.append((ap, bp))
        try:
            self.a_index[ap].append(bp)
        except:
            self.a_index[ap] = [bp]
        try:
            self.b_index[bp].append(ap)
        except:
            self.b_index[bp] = [ap]

    def __collide(self, A, B):
        """Actual, common ref, collision test"""

        bstart = ap = bp = 0

        while(ap < len(A)):
            first = None
            a = A[ap]
            b = B[bp]

            while(a.Locus().stop >= b.Locus().start):
                if(b.Locus().stop >= a.Locus().start):
                    self.__overlap(a, b)
                    if(first == None):
                        first = bp
                bp += 1
                if(bp >= len(B)):
                    break
                b = B[bp]

            ap += 1
            if(first != None):
                bstart = first
            bp = bstart

    def __getitem__(self, i):
        """Return the ith overlap as a tuple of (A index, B index)"""
        return self.__overlaps[i]

    def __getslice__(self, i, j):
        return self.__overlaps[i:j]

    def __len__(self):
        return len(self.__overlaps)

    def OverlapA(self, i):
        """Return a list of elements from B overlapping i, where
        i is an element of A."""
        try:
            return self.a_index[i]
        except:
            return []

    def OverlapB(self, i):
        """Return a list of elements from A overlapping i, where
        i is an element of B."""
        try:
            return self.b_index[i]
        except:
            return []

def LocusDiffDisjoint(lhs, rhs, presorted = False):
    """Return the maximum set of loci spanned by lhs and not by rhs as
    for LocusDiff but for the special case where rhs is known to be
    disjoint (no overlapping pairs).  If presorted=True then rhs is
    known to be sorted by start."""
    if(not presorted):
        rhs = sorted(rhs, key = lambda x: x.start)
    loci = []
    locus = lhs.Locus()
    if(locus.start < rhs[0].Locus().start):
        loci.append(Locus.fromPrototype(locus, stop = rhs[0].start-1))
    for n in range(len(rhs)-1):
        start = rhs[n].stop+1
        stop = rhs[n+1].start-1
        if(start <= stop):
            loci.append(Locus.fromPrototype(locus, start=start, stop=stop))
    if(locus.stop > rhs[-1].Locus().stop):
        loci.append(Locus.fromPrototype(locus, start=rhs[-1].stop+1))

    return loci
        
def LocusDiff(lhs, rhs):
    """Return the maximum set of loci spanned by lhs and not by rhs where:
    lhs.Locus() is defined.
    rhs is an iterable of objects with Locus() defined.
    (strand is not considered for this calculation)

    The returned value is a list with its disjoint elements in cmpstart
    order.
    """

    # Implementing in terms of bitmap a la ClassifyBases.
    #   The bitmask implementation is definitely simpler than scanning a sorted
    #      list (fewer corner cases, fewer steps, ...)
    #   It definitely requires more memory (but potentially fewer intermediate
    #       data structures).
    #   I'm not sure which method is faster...

    # Note: absent a bitmask class in python, we pay 1 byte per bit.
    locus = lhs.Locus()
    bases = [0]*len(locus)
    for i in rhs:
        r = i.Locus()
        if((r.ref != locus.ref) or
           (r.start > locus.stop) or
           (r.stop < locus.start)):
            continue
        for j in range(max(locus.start, r.start) - locus.start,
                       min(locus.stop, r.stop) - locus.start + 1):
            bases[j] |= 1

    p0 = 0
    retval = []
    while(p0 < len(bases)):
        while((p0 < len(bases)) and (bases[p0] != 0)):
            p0 += 1
        if(p0 >= len(bases)):
            break
        p1 = p0
        while((p1 < len(bases)) and (bases[p1] == 0)):
            p1 += 1
        retval.append(Locus(ref = locus.ref,
                            start = p0 + locus.start,
                            stop = p1-1 + locus.start,
                            strand = locus.strand,
                            genome = locus.genome))

        p0 = p1

    return retval


if(__name__ == "__main__"):
    class dummy:
        def __init__(self, start, stop):
            self.ref = None
            self.start = start
            self.stop = stop
        def __repr__(self):
            return "(%d, %d)" % (self.start, self.stop)
        def Locus(self):
            return self
    
    A = list(map(lambda x: dummy(x[0],x[1]),
                 ([10,20], [30,40], [50,90], [52,100], [55,70])))
    B = list(map(lambda x: dummy(x[0],x[1]),
                 ([5,25], [20,45], [55,70], [60, 120])))

    print(A)
    print(B)
    overlaps = Collisions(A,B)
    print(list(map(lambda x: (A[x[0]], B[x[1]]), overlaps)))

    overlaps = RefCollisions(A,B)
    print(list(map(lambda x: x, overlaps)))

