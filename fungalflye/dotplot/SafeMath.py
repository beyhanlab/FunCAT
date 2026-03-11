#!/usr/bin/env python3
"""
Math functions that can handle None values.
Primarily intended for math on microarray log
ratios with missing values.
"""

from math import log
log2 = log(2)
def safelog(ratio,base=2):
    if(base == 2):
        b = log2
    else:
        b = log(base)
    try:
        return log(float(ratio))/b
    except:
        return None
def safeadd(a,b):
    try:
        return a+b
    except:
        return None
def safesub(a,b):
    try:
        return a-b
    except:
        return None
def safesum(x):
    s = None
    for i in x:
        if(i is not None):
            if(s is None):
                s = i
            else:
                s += i
    return s
def safemean(x):
    s = 0.
    c = 0
    for i in x:
        if(i is not None):
            s += i
            c += 1
    if(c == 0):
        return None
    return s/c

if(__name__ == "__main__"):
    print("Hello, world")
