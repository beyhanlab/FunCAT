#!/usr/bin/env python3
# Generic utility functions

def transpose(a):
    """Given rectangular list of lists a, return its transpose as a
    list of lists."""
    return [[j[i] for j in a] for i in range(len(a[0]))]

def gapped_fp(fp, skips):
    """Yield all lines in a filestream except for line numbers
    in skips (counting from zero)."""
    for (i,line) in enumerate(fp):
        if(i in skips):
            continue
        yield line

def bin_to_text(fp):
    """Simple wrapper to autodecode binary streams to text"""
    for line in fp:
        yield line.decode()
        
def hexdump(data, offset, lines):
    """Return data[offset:offset+lines*16] formatted like hexdump -C"""
    for i in range(lines):
        curoff = offset+16*i
        return "%08x  %s  %s  |%s|" % (
            curoff,
            " ".join("%02x" % ord(data[curoff+j]) for j in range(8)),
            " ".join("%02x" % ord(data[curoff+j]) for j in range(8,16)),
            # TODO: replace dots with characters for "printable" characters
            "................")

def array_index(a,i):
    """Return index to first occurence of i in sequence a.
    Intent is to provide the equivalent of list.index for numpy arrays."""
    for (n,j) in enumerate(a):
        if(j == i):
            return n
    raise ValueError("%s is not in array" % i)
    
# N.B. Removed function base2string(i) was equivalent to builtin bin(i)[2:]

# Naive median function for when R is not available or not cooperating.
# Note that we could probably get away with a partial sort.
def median(s, remove_none = True):
    """Return the median value of iterable s, assuming that
    there is a meaningful rank-ordering of s's contents."""
    if(remove_none):
        s = [i for i in s if(i is not None)]
    t = sorted(s)
    if(len(t) == 0):
        return None
    p = int(len(t)/2.0)
    if(len(t) % 2 == 0):
        return float(t[p-1]+t[p])/2.0
    else:
        return t[p]

# MySQL short-hand
def get_last_insert_id(cursor):
    import MySQLdb
    try:
        cursor.execute("SELECT LAST_INSERT_ID()")
        return cursor.fetchone()[0]
    except MySQLdb.Error as e:
        sys.stderr.write("Error getting insert_id!\n%s\n" % e)

# Decorator to add a "prototype" constructor to a class
#  Usage:
#   from MsvUtil import prototypeable
#   @prototypeable
#   class MyClass:
#      def __init__(self, ...):
#         ...
#   a = MyClass(...)
#   b = MyClass.fromPrototype(a, ...)
#
#  Note that this currently fails for classes with member data that is
#    not taken directly from __init__ parameters

def fromPrototype(cls, prototype, **kw):
    for (key,val) in prototype.__dict__.items():
        if(not (key in kw)):
            kw[key] = val

    return cls(**kw)

def prototypeable(cls):
    cls.fromPrototype = classmethod(fromPrototype)
    return cls

# Based on http://docs.python.org/dev/library/itertools.html
# (made docstring more explicit and backported from Python3 to Python2.6
#  version of itertools)
# Worst case behavior is to buffer all of one result while the other is
# being processed, which is not that different from the memory consumption
# of generating both lists explicitly.
def partition(pred, iterable):
    """Return (false, true) pair of iterators partitioning iterable
    on pred(icate)."""
    from itertools import tee, ifilter, ifilterfalse
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = tee(iterable)
    return ifilterfalse(pred, t1), ifilter(pred, t2)

class ParenIndex:
    """An bi-directional index of matching left/right parens on a string."""
    def __init__(self, s, left = "{", right = "}"):
        self.pairs = []
        stack = []
        for (i,c) in enumerate(s):
            if(c == left):
                stack.append(i)
            elif(c == right):
                self.pairs.append((stack[-1],i))
                stack = stack[:-1]
        assert(len(stack) == 0)
        self.left_index = dict(self.pairs)
        self.right_index = revdict(self.left_index)

    def find_after(self, i):
        """Return first pair with left paren at or after i.
        Raise IndexError if there is no such pair."""
        for p in sorted(self.pairs):
            if(p[0] >= i):
                return p
        raise IndexError

    def find_enclosing(self, i):
        """Return pairs spanning index i, sorted parent to child."""
        r = []
        for p in self.pairs:
            if(p[0] <= i <= p[1]):
                r.append(p)
        r.sort()
        return r

# N.B.: cmp_index and cmp_f can usually be replaced with the key argument
#       of the standard sort function.

class cmp_index:
    """Class for comparing indexable containers on one or more indices."""
    # TODO: extend to arbitrary number of indices
    def __init__(self, i):
        # Iterable check from
        #  http://mail.python.org/pipermail/python-list/2006-July/394493.html
        try:
            it = iter(i)
            self.i = i
        except TypeError:
            self.i = (i,)
            
    def __call__(self, lhs, rhs):
        for indice in self.i:
            t = lhs[indice].__cmp__(rhs[indice])
            if(t):
                return t
        return 0

class cmp_f:
    """Class for comparing objects on one or more observables."""
    def __init__(self, i):
        """Init from a callable or iterable of callables."""
        # Iterable check from
        #  http://mail.python.org/pipermail/python-list/2006-July/394493.html
        try:
            it = iter(i)
            self.i = i
        except TypeError:
            self.i = (i,)
            
    def __call__(self, lhs, rhs):
        for f in self.i:
            t = cmp(f(lhs),f(rhs))
            if(t):
                return t
        return 0

def f_args(f):
    """Return the names of the arguments of function f as two tuples: a
    simple tuple of required arguments and a tuple of (name, default)
    for optional arguments.
    """
    args = f.__code__.co_varnames
    defaults = a.func_defaults
    n = len(args)-len(defaults)
    return tuple(args[:n]),tuple(zip(args[n:],defaults))

def f_defaults(f):
    """Return a dictionary of the default arguments of function f."""
    return dict(get_args(f)[1])
    
def revdict(d):
    """Return a reverse mapping for the 1-to-1 mapping represented by
    d (which supports the dict.items() interface).
    """
    return dict((b,a) for (a,b) in d.items())

def multirevdict(d):
    """Return a reverse mapping for the many-to-1 mapping represented by
    d (which supports the dict.items() interface).
    """
    retval = {}
    for (key, val) in d.items():
        try:
            retval[val].append(key)
        except KeyError:
            retval[val] = [key]

    return retval

def hdict(s, f = lambda x: x):
    """Return a "histogram" dictionary for a sequence or
    function of a sequence.

    Note that this is equivalent to indexing s on f(s[i])"""
    retval = {}
    for (key, val) in ((f(j), j) for j in s):
        try:
            retval[key].append(val)
        except KeyError:
            retval[key] = [val]

    return retval

def argmax(s, f = lambda x: x):
    """Return the element of sequence s corresponding to the first
    maximum value of f(s[i]).
    Return None for empty sequence (note that sequences with None
    elements have ambiguous return values)."""
    retval = None
    best = None
    for i in s:
        val = f(i)
        if((best is None) or (val > best)):
            best = val
            retval = i

    return retval

def argmin(s, f = lambda x: x):
    """Return the element of sequence s corresponding to the first
    minimum value of f(s[i]).
    Return None for empty sequence (note that sequences with None
    elements have ambiguous return values)."""
    retval = None
    best = None
    for i in s:
        val = f(i)
        if((best is None) or (val < best)):
            best = val
            retval = i

    return retval

# CSV wrapper class

class Table:
    def __init__(self, header, rows):
        self.header = header
        self.rows = rows
        self.name2col = dict((i,n) for (n,i) in enumerate(self.header))
        
    def __getitem__(self, i):
        if(isinstance(i, slice)):
            return self.getslice(i.start, i.stop)
        if(isinstance(i, int)):
            if(0 <= i <  len(self.rows)):
                return TableRow(self, i)
            elif(-len(self.rows) <= i < 0):
                return TableRow(self, len(self)+i)
            else:
                raise IndexError
        else:
            try:
                j = self.name2col[i]
                return [k[j] for k in self.rows]
            except KeyError:
                # N.B. Changing exception type for backwards compatibility
                #      with old implementation.
                raise IndexError

    def getslice(self, i, j):
        if(i is None):
            i = 0
        if(j > len(self)):
            j = len(self)
        assert((type(i) == type(j) == int))
        assert(0 <= i < len(self))
        assert(0 <= j <= len(self))
        # TODO: support negative indices
        return Table(header = self.header,
                     rows = self.rows[i:j])

    def keys(self):
        return self.name2col.keys()
    
    def __len__(self):
        return len(self.rows)

    def pad(self, null = ""):
        """Pad ragged rows to match length of header with null value.
        Raise ValueError if any row is longer than the header."""
        for i in self.rows:
            n = len(self.header) - len(i)
            if(n < 0):
                raise ValueError
            if(n > 0):
                i += [null]*n

    @classmethod
    def fromCsv(cls, fin, skip = 0, comment = None):
        """Generate a Table from a CSV file.
        fin is a filename or an object
            supporting the file interface.
        """
        from csv import reader
        if(type(fin) == str):
            fp = reader(open(fin,"rt"))
        else:
            fp = reader(fin)

        for i in range(skip):
            dummy = next(fp)
            
        return cls(header = next(fp),
                   rows = [i for i in fp
                           if((comment is None) or
                              (len(i) == 0) or
                              (not i[0].startswith(comment)))])

    def writeCsv(self, fout):
        """Write the table to fout in Excel-compatible CSV format."""
        from csv import writer
        out = writer(fout)
        out.writerow(self.header)
        for row in self:
            out.writerow(list(row))

    @classmethod
    def fromExcelTdt(cls, fin, comment = None, skip = 0):
        """Generate a Table from a TDT file.
        fin is a filename or an object
            supporting the file interface.
        Uses Excel rules for things like quoting.
        Skip lines starting with comment, if given.
        """
        from csv import reader, excel_tab
        if(type(fin) == str):
            fp = reader(open(fin,"rt"), dialect = excel_tab)
        else:
            fp = reader(fin, dialect = excel_tab)

        for i in range(skip):
            dummy = next(fp)

        if(comment is None):
            return cls(header = next(fp), rows = [i for i in fp])

        else:
            header = next(fp)
            while(header[0].startswith(comment)):
                header = next(fp)

            return cls(header = header,
                       rows = [i for i in fp if(not i[0].startswith(comment))])
        
    @classmethod
    def fromTdt(cls, fin, skip = 0, comment = None):
        """Generate a Table from a TDT file.
        fin is a filename or an object
            supporting the file interface.
        This is a naive parse of \t delim with \r\n or \n EOLs
        and does not support field quoting.
        """
        
        if(type(fin) == str):
            fp = open(fin,"rt")
        else:
            fp = fin

        for i in range(skip):
            dummy = next(fp)
            
        return cls(header = next(fp).rstrip("\r\n").split("\t"),
                   rows = [line.rstrip("\r\n").split("\t") for line in fp
                           if((comment is None) or
                              (not line.startswith(comment)))])

    @classmethod
    def fromMacTdt(cls, fin):
        """Generate a Table from a TDT file with \r EOLs.
        fin is a filename or an object
            supporting the file interface.
        This is a naive parse of \t delim
        and does not support field quoting.
        """
        
        if(type(fin) == str):
            fp = open(fin,"rt")
        else:
            fp = fin

        lines = fp.read().decode().split("\r")

        return cls(header = lines[0].split("\t"),
                   rows = [line.split("\t") for line in lines[1:]])

    def writeTdt(self, fout):
        """Write the table to fout in Excel-compatible Tdt format."""
        from csv import writer, excel_tab
        out = writer(fout, dialect = excel_tab)
        out.writerow(self.header)
        for row in self:
            out.writerow(list(row))

    @classmethod
    def fromXls(cls, fin):
        """Generate a Table from an XLS (Excel binary) file.
        fin is a filename.
        """

        from subprocess import Popen, PIPE
        p = Popen(("xls2csv",fin), stdout = PIPE)
        retval = cls.fromCsv(p.stdout)
        # xls2csv appends a line-feed character at the end of each
        #   page.  Currently, we assume only one page and kill the
        #   final linefeed.  Alternatively, we could split on linefeeds
        #   and generate a separate Table instance for each page.
        assert(retval.rows[-1][0] == '\x0c')
        retval.rows = retval.rows[:-1]
        return retval

    def varcol(self):
        """Print summary of unique elements in each table column."""
        print("Total rows:", len(self))
        for (n,i) in enumerate(self.header):
            print("%s: %d" % (i, len(set(j[n] for j in self))))

class TableRow:
    def __init__(self, table, row):
        self.table = table
        self.row = row

    def dict(self):
        """Return a dict of {colname => value} for this row.

        N.B.: TableRow implements __getitem__ and keys(), which
              makes it a python mapping.  Therefore, it is usually
              preferable to treat a TableRow as a dict; 
              e.g.:
                x = row["foo"]
                "{foo}.{bar}".format(**row)
              
              The dict method is provided for cases where an actual
              independent dict instance is required.  E.g., when
              using this row as the starting point for a dict to
              be augmented with additional keys (c.f. RnaSeqPipeline.ipynb).
        """
        return dict(zip(self.table.header,
                        self.table.rows[self.row]))

    def keys(self):
        return self.table.keys()
    
    def __str__(self):
        return "\n".join("%s:\t\t%s" % i
                         for i in zip(self.table.header,
                                      self.table.rows[self.row]))
    def __repr__(self):
        return str(self)
    
    def __getitem__(self, i):
        if(isinstance(i, slice)):
            return self.getslice(i.start, i.stop)
        if(type(i) == int):
            return self.table.rows[self.row][i]
        else:
            try:
                return self.table.rows[self.row][self.table.name2col[i]]
            except KeyError:
                # N.B. Changing exception type for backwards compatibility
                #      with old implementation.
                raise IndexError

    def getslice(self, i, j):
        """Return row values for columns i through j.
        i and j can be given (independently) as indices or
        column names.  Indices are interpreted as for python
        lists (first index is 0, i is the index of the first
        returned element, j is one past the index of the last
        returned element).  Column names are interpreted as
        for R data.frames (i is the column of the first returned
        element, j is the column of the last returned element).
        """
        if(type(i) == int):
            ii = i
        elif(i is None):
            ii = None
        else:
            try:
                ii = self.table.name2col[i]
            except KeyError:
                raise IndexError

        if(type(j) == int):
            jj = j
        elif(j is None):
            jj = None
        else:
            try:
                jj = self.table.name2col[j]+1
            except KeyError:
                raise IndexError

        return self.table.rows[self.row][ii:jj]

    def __len__(self):
        try:
            return len(self.table.rows[self.row])
        except IndexError:
            print(self.row)
            raise

class IterTable:
    """Partial implementation of Table for streaming.
    TODO: Table should inherit from this class or a common abstract base
    """
    def __init__(self, header, rowiter):
        self.header = header
        self.rowiter = rowiter
        self.name2col = dict((i,n) for (n,i) in enumerate(self.header))

    def __iter__(self):
        for row in self.rowiter:
            yield IterTableRow(self, row)

class IterTableRow:
    def __init__(self, table, fields):
        self.table = table
        self.fields = fields

    def dict(self):
        """Return a dict of {colname => value} for this row.

        N.B.: TableRow implements __getitem__ and keys(), which
              makes it a python mapping.  Therefore, it is usually
              preferable to treat a TableRow as a dict; 
              e.g.:
                x = row["foo"]
                "{foo}.{bar}".format(**row)
              
              The dict method is provided for cases where an actual
              independent dict instance is required.  E.g., when
              using this row as the starting point for a dict to
              be augmented with additional keys (c.f. RnaSeqPipeline.ipynb).
        """
        return dict(zip(self.table.header,
                        self.fields))

    def keys(self):
        return self.table.keys()
    
    def __str__(self):
        return "\n".join("%s:\t\t%s" % i
                         for i in zip(self.table.header,
                                      self.fields))
    def __repr__(self):
        return str(self)
    
    def __getitem__(self, i):
        if(isinstance(i, slice)):
            return self.getslice(i.start, i.stop)
        elif(isinstance(i, int)):
            return self.fields[i]
        else:
            try:
                return self.fields[self.table.name2col[i]]
            except KeyError:
                # N.B. Changing exception type for backwards compatibility
                #      with old implementation.
                raise IndexError

    def getslice(self, i, j):
        """Return row values for columns i through j.
        i and j can be given (independently) as indices or
        column names.  Indices are interpreted as for python
        lists (first index is 0, i is the index of the first
        returned element, j is one past the index of the last
        returned element).  Column names are interpreted as
        for R data.frames (i is the column of the first returned
        element, j is the column of the last returned element).
        """
        if(type(i) == int):
            ii = i
        elif(i is None):
            ii = None
        else:
            try:
                ii = self.table.name2col[i]
            except KeyError:
                raise IndexError

        if(type(j) == int):
            jj = j
        elif(j is None):
            jj = None
        else:
            try:
                jj = self.table.name2col[j]+1
            except KeyError:
                raise IndexError

        return self.fields[ii:jj]

    def __len__(self):
        try:
            return len(self.fields)
        except IndexError:
            print(self.fields)
            raise
        
class FacetedIndex:
    """Experimental class for indexing a multidimensional associative array
    in arbitrary order.  Assumes that keys are unique across all dimensions.
    C.f. NickMerge.py for an example of using this class."""
    def __init__(self, dimensions):
        if(isinstance(dimensions,int)):
            self.dimensions = dimensions
            self.dimnames = ["V%d" % i for i in range(self.dimensions)]
        else:
            self.dimensions = len(dimensions)
            self.dimnames = tuple(dimensions)
        self.data = {}
        self.key2dim = {}
    def put(self, path, val):
        assert(len(path) == self.dimensions)
        d = self.data
        for (n,i) in enumerate(path[:-1]):
            try:
                d = d[i]
            except KeyError:
                d[i] = {}
                d = d[i]
                self.key2dim[i] = n
        d[path[-1]] = val
        self.key2dim[path[-1]] = self.dimensions-1
    def get(self, path):
        d = self.data
        for i in path:
            d = d[i]
        return d
    def __getitem__(self, i):
        s = [None]*self.dimensions
        s[self.key2dim[i]] = i
        return FacetProxy(self, s)
        
class FacetProxy:
    def __init__(self, parent, s):
        self.parent = parent
        self.s = s[:] 
    def __getitem__(self,i):
        d = self.parent.key2dim[i]
        assert(self.s[d] is None)
        s = self.s[:d]+[i]+self.s[d+1:]
        if(None in s):
            return FacetProxy(self.parent,s)
        else:
            return self.parent.get(s)

# Unit test
if(__name__ == "__main__"):
    # N.B. for proper table testing, we need some test files in git as well
    print("Hello, world")


