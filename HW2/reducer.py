#!/usr/bin/env python

import sys

# Initialize
current_key = None
current_count = 0

# input comes from STDIN
for line in sys.stdin:
    # remove leading and trailing white space
    line = line.strip()
    
    # split the input we get from mapper.py
    # into key and count
    key, count = line.split('\t')
    
    # convert count to type integer
    count = int(count)
    
    # this IF-switch only works because Hadoop sorts map output
    # by key before it is passed to the reducer
    # if the new key is the same as the current key,
    # increment the current count by 1
    # or more generally, by count
    if current_key == key:
        current_count += count
    # if the new key is different from the current key...
    else:
        # if the current key is NOT None...
        if current_key:
            # write result to STDOUT
            print '%s,%s' % (current_key,current_count)
        # if the new key is different from the current key
        # reset the current count to 1
        # or more generally, to count
        # and set the current key to the new key
        current_count = count
        current_key = key

# do not forget to output the last result if needed!
if current_key == key:
    print '%s,%s' % (current_key,current_count)