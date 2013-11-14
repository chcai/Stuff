#!/usr/bin/env python

# import 'packages'
import sys
import math

# input comes from STDIN (standard input)
for line in sys.stdin:
    # remove leading and trailing whitespace
    line = line.strip()
    
    # split the line into x and y
    entries = line.split('\t')
    
    # key will be 'x_lo,x_hi,y_lo,y_hi'
    # calculate x_lo and y_lo
    # we want bins for x and y with length .1
    # so we want only 1 decimal place
    x_lo=math.floor(float(entries[0])*10)/10
    y_lo=math.floor(float(entries[1])*10)/10
    
    # spit out (key, value)
    # where key and value are tab seperated
    # and within key, nothing is tab separated
    # since we are counting, the value is a count, 1
    print '%s,%s,%s,%s\t%s' % (x_lo,x_lo+.1,y_lo,y_lo+.1,1)