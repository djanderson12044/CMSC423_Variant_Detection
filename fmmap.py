# Any and all necessary documentation is in the README (if there is any)

# TODO Important notes
# 1. the great python cprofile profiler. It will profile your call stack to let you know, empirically,
#  how time is being spent in the functions within your program. It might be exactly where you expect, but it might be
#  somewhere else entirely.

import math
import gzip
import sys
from Bio import SeqIO
import pickle  # Writes and reads binary foramt


def main():
    # Just some basic input checking
    if len(sys.argv) > 1:
        if sys.argv[1] == "index":
            if len(sys.argv) != 4:
                print("Error: Wrong Number of Arguments for Index")
            else:
                index("data/" + sys.argv[2], "data/" + sys.argv[3])
        elif sys.argv[1] == "align":
            if len(sys.argv) != 5:
                print("Error: Wrong Number of Arguments for Align")
            else:
                align("data/" + sys.argv[2], "data/" + sys.argv[3], "data/" + sys.argv[4])
    else:
        print("Error: Wrong Number of arguments")


# ref: (.fa) file that needs to be indexed.
# ref_index: (.P)output file to store all the index information
def index(ref, ref_index):
    sa = []  # Suffix Array

    lColumn = []  # Last column of BWM
    fColumn = []  # First column of BWM
    occ = {}  # Rankings, dictionary of arrays

    # Checks to see if the reference file exists in data
    try:
        record = SeqIO.read(ref, "fasta")  # biopython read of the ref file
    except FileNotFoundError:
        print("Reference file not found")
        return

    text = list(record.seq)  # made it so that I could add the $ terminal value
    text.append('$')

    return_index = {}  # What we are writing to ref_index in binary format using pickle

    # Populate Suffix Array
    for x in range(0, len(text)):
        sa.append(x)

    quicksort(sa, text, 0, len(text) - 1)

    # Creates last column from suffix array
    # And creates occ from last column -------------------------
    for n in range(len(sa)):
        c = text[sa[n] - 1] if sa[n] > 0 else '$'

        lColumn.append(c)

        # Adding to occ
        if c in occ:
            pass
        else:
            if c != '$':
                occ[c] = []
                temp = c

                for k in occ.keys():
                    if k != temp:
                        temp = k
                        break

                for _ in range(len(occ[temp])):
                    occ[c].append(0)

        for x in occ.keys():
            lastElement = occ[x][-1] if len(occ[x]) > 0 else 0
            if x == c:
                occ[x].append(lastElement + 1)
            else:
                occ[x].append(lastElement)
        # OCC END -------------------------------------------

    # Creates first column
    fColumn.append('$')
    for n in sorted(occ.keys()):
        for x in range(occ[n][-1]):
            fColumn.append(n)

    # Adding to dictionary which is converted to binary at the end of the program
    return_index["name"] = record.id
    return_index["length"] = len(text)
    return_index["sa"] = sa
    return_index["occ"] = occ
    return_index["lColumn"] = lColumn
    return_index["fColumn"] = fColumn
    return_index["ref"] = record.seq  # The entire reference string

    with open(ref_index + ".P", "w+b") as outfile:
        pickle.dump(return_index, outfile)
        outfile.close()

    print(record.seq)


# ref_index: (.?) File from the index definition
# reads: (.fa)
# alignments: (.sam)
# Mapping Algorithm
def align(ref_index, reads, aligns):  # todo Finish implementing align method
    ninf = float("-inf")
    #seed_skip = lambda l: math.floor(l / 5.0) # TODO Real Seed Skip Function
    seed_skip = lambda l: math.floor(l / 2.0)
    gap = 5

    # with gzip.open(reads + '.fa.gz', 'rt') as rfile: # TODO Real file is gziped
    for read in SeqIO.parse(reads, "fasta"):
        alignments = []
        read_len = len(read.seq)
        best_score = ninf
        seed_pos = 0
        skip = seed_skip(read_len)
        for seed_start in range(0, read_len, skip):
            seed_end = min(read_len, seed_start + skip)
            interval, match_len = get_interval(read.seq[seed_start:seed_end], ref_index)


        # given all the places where the seed matches, look for an alignment around it
        # the ref_positions member of `bwt_index` will return positions on the reference
        # string corresponding to the *beginning of the read*, assuming there are no gaps
        # in the alignment before the start of the seed (handling that is why we do fitting
        # alignment below).
        for ref_pos in ref_positions(interval, seed_end, match_len):
            pass

        '''
            # Perform a "fitting" alignment of the query (seq) into the reference (ref)
            # the returned alignment object contains the score, the alignment and the
            # implied position where the query (seq) begins under the alignment.
            # To perform the fitting alignment, you should "slice out" a region of the
            # reference around the implied start position (ref_pos) with a bit of padding
            # (e.g. gap bases) before the first base where the read would start and after
            # the last base where the read would end.  This will ensure that the fitting_alignment
            # procedure can find the optimal position for the query within the reference
            # window that contains it, even if there are insertions or deletions in the read.
            alignment = fitting_alignment(seq, ref, ref_pos, gap)
            if alignment.score > best_score:
                best_score = alignment.score
                alignments = [alignment]
            elif alignment.score == best_score:
                alignments.append(alignment)
        for a in alignments:
            write_to_sam(output_file, a)
        '''


# Align Helper Methods-----------------------------------------

def get_interval(pattern, ref_index):
    r_index = {}
    with open(ref_index + ".P", "r+b") as infile:
        r_index = pickle.load(infile)

    # x is start of interval, y is end of interval
    x, y = bbwm(r_index["occ"], r_index["lColumn"], pattern)

    print(pattern)
    print(x, y)
    print("--------------------------")

    for n in range(x, y + 1):
        print(n, r_index["ref"][r_index["sa"][n]:])

    print()

    if x < 0 or y < 0:
        return (-1, -1), -1
    else:
        return (x, y), len(pattern)


def ref_positions(interval, seed_end, match_len): # TODO Need to implement
    return 1, 2


# Helper Methods----------------------------------------

# Using quick sort algorithm
# Inspiration: https://medium.com/human-in-a-machine-world/quicksort-the-best-sorting-algorithm-6ab461b5a9d0
# arr: unordered SA
# string: what we are creating the SA frm
# left: At the beginning is 0
# right: At the beginning is len(string)
def quicksort(arr, string, l, r):
    if l >= r:
        return

    pivot = r
    count = l

    for i in range(l, r + 1):
        if string[arr[i]:] <= string[arr[pivot]:]:
            temp = arr[count]
            arr[count] = arr[i]
            arr[i] = temp
            count += 1

    quicksort(arr, string, l, count - 2)
    quicksort(arr, string, count, r)


# Back word searching using the the fm-index
# Returns the range of the SA at which the pattern was matched to the reference
def bbwm(occ, lc, p):
    top = 0  # Start of the range
    bottom = len(lc) - 1  # End of the range

    pattern = list(p)  # Seed that we are matching against the reference

    # Initializes the count table, which we can query to get the number of characters
    #  appearing before the character set that we want (in the first column)
    curr = 1;
    count_table = {}
    count_table['$'] = 0  # Not stored in occ so we have to add it
    for x in sorted(occ.keys()):
        count_table[x] = curr
        curr += occ[x][-1]

    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern.pop()

            if symbol in lc[top: bottom + 1]:
                top = count_table[symbol] + (occ[symbol][top - 1] if top != 0 else 0)
                bottom = count_table[symbol] + occ[symbol][bottom] - 1
            else:
                return -1, -1
        else:
            return top, bottom
    return -1, -1


# Calls the main() definition
main()

# Command Line parameters
# "index" "test_ref.fa" "ref_index"
# "align" "ref_index" "reads_simple_copy.fa" "alignments"
