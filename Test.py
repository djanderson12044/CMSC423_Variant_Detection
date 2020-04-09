# RPS 10
import sys
import os
import copy


# n: symbol
# i: stop index, up to but not including
# lc: last column
# Gets the rank of n in the slice from lc[0:i]
def count(n, i, lc):
    count = 0
    for x in range(0, i):
        if n == lc[x]:
            count += 1

    return count


# fo: firstoccurence (array of indices of the first characters appearing
# in first column. Works in tandem with ot.)
# ot: order tracker (sorted list of all the chars in the text)
# lc: last column
# p: pattern
# c: count
def bbwm(fo, ot, lc, p):
    top = 0
    bottom = len(lc) - 1
    while top <= bottom:
        if len(p) > 0:
            symbol = p[-1]
            # removing last letter from pattern
            temp = list(p)
            temp.pop()
            p = "".join(temp)

            if symbol in lc[top: bottom + 1]:
                top = fo[ot.index(symbol)] + count(symbol, top, lc)
                bottom = fo[ot.index(symbol)] + count(symbol, bottom + 1, lc) - 1
            else:
                return 0
        else:
            return top, bottom


# Above are helper methods-------------

# Openes data set and output file
file = open(sys.argv[1], "r")

lastcolumn = file.readline()
lastcolumn = lastcolumn.rstrip()  # Strips new line character

patterns = file.readline()
patterns = patterns.rstrip()
# really only thing that matters
patternlist = patterns.split(" ")

# Works in tandem with ordertracker
firstoccurrence = list()

# How many of each character
chartracker = dict()
# Sorted array of what characters are in text
ordertracker = list()

for n in lastcolumn:
    if n in chartracker:
        chartracker[n] += 1
    else:
        chartracker[n] = 1
        ordertracker.append(n)

ordertracker.sort()

# Making first occurrence array
curr = 0;
for x in ordertracker:
    firstoccurrence.append(curr)
    curr += chartracker[x]

result = list()

for n in range(0, len(patternlist)):
    x = bbwm(firstoccurrence, ordertracker, lastcolumn, patternlist[n])
    print(x)

file.close()
