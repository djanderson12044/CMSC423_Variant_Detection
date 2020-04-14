# Any and all necessary documentation is in the README (if there is any)

import math
import gzip
import sys
from Bio import SeqIO
import pickle  # Writes and reads binary format
import pysam # Don't actually use it

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


# ref_index: (.?) File from the index definition
# reads: (.fa)
# alignments: (.sam)
# Mapping Algorithm
def align(ref_index, reads, aligns):  # todo Finish implementing align method
    ninf = float("-inf")
    seed_skip = lambda l: math.floor(l / 5.0) # TODO Real Seed Skip Function
    gap = -2  # Cost for insertion/

    # File we are writing to inplace of the sam formatted file
    output = open(aligns + ".txt", "w")

    r_index = {}
    with open(ref_index + ".P", "r+b") as infile:
        r_index = pickle.load(infile)

    # Include and indent if trying to read from a gzipped file
    #with gzip.open(reads + '.fa.gz', 'rt') as rfile: # TODO Real file is gziped
    for read in SeqIO.parse(reads, "fasta"):
        alignments = []
        read_len = len(read.seq)
        best_score = ninf
        seed_pos = 0
        skip = seed_skip(read_len)
        q_len = len(read.seq)

        for seed_start in range(0, read_len, skip): # TODO need to put this back
            seed_end = min(read_len, seed_start + skip) # offset of the read where the seed ends.

            interval, match_len = get_interval(read.seq[seed_start:seed_end], ref_index, r_index["occ"], r_index["lColumn"])

            for ref_pos in ref_positions(interval, seed_end, match_len, q_len, r_index["sa"]):
                alignment = fitting_alignment(read.seq, r_index["ref"], ref_pos, gap)
                if alignment["score"] > best_score: #
                    best_score = alignment["score"]
                    alignments = [alignment]
                elif alignment["score"] == best_score:
                    alignments.append(alignment)
            for a in alignments:
                output.write("read_name:" + read.id + "\r\n")
                output.write("Query: " + "".join(read.seq) + "\r\n")
                output.write("Cigar:" + a["cigar"] + "\r\n")
                output.write("Score:" + str(a["score"]) + "\r\n")
                output.write("Position in Ref: " + str(a["position"]) + "\r\n")
                output.write("-------------------------------\r\n")
            #    write_to_sam(output_file, a)
    output.close


# Align Helper Methods-----------------------------------------
def get_interval(pattern, ref_index, occ, lcolumn):
    # x is start of interval, y is end of interval
    x, y = bbwm(occ, lcolumn, pattern)

    if x < 0 or y < 0:
        return (-1, -1), -1
    else:
        return (x, y), len(pattern)


# interval: interval of sa indexes that match exactly with the pattern
# sed_end: offset of the read where the seed ended
# match_len: length of the match
def ref_positions(interval, seed_end, match_len, q_len, sa):
    positions = []

    if match_len > 0:
        for p in range(interval[0], interval[1] + 1):
            temp = sa[p] - (seed_end - match_len)
            if temp - 5 >= 0 and temp + q_len + 5 < len(sa) - 1:
                positions.append(temp)

    return positions


# The returned alignment object contains the score, the alignment and the
#  implied position where the query (seq) begins under the alignment.
def fitting_alignment(seq, ref, ref_pos, gap):
    ref_slice = ref[ref_pos - 5: ref_pos + len(seq) + 5]
    match = 0
    mismatch = -2
    insertion = gap
    deletion = gap

    result = {}

    opt = [[]]
    # -5: terminal
    # 1: left
    # 0: diagonal
    # -1: down
    d_opt = [[]]  # Stores direction of the opt

    max_score = float("-inf")
    max_opt = (-1, -1)

    # Building opt table ------------------------------
    for x in range(0, len(ref_slice) + 1):
        opt[0].append(0)
        d_opt[0].append(-5)

    for j in range(1, len(seq) + 1):
        opt.append([j * gap])
        d_opt.append([-1])
        for i in range(1, len(ref_slice) + 1):
            diag = (match if ref_slice[i - 1] == seq[j - 1] else mismatch) + opt[j - 1][i - 1]
            try:
                left = insertion + opt[j][i - 1]
            except IndexError:
                print("Here is the left error")
            try:
                down = deletion + opt[j - 1][i]
            except IndexError:
                print("Here is the down error")

            temp = max(diag, left, down)

            # Initializing d_opt table
            if diag == temp:
                d_opt[j].append(0)
            elif left == temp:
                d_opt[j].append(1)
            elif down == temp:
                d_opt[j].append(-1)
            else:
                print("OH FRICK, something went wrong")

            if j == len(seq):
                if temp >= max_score:
                    max_score = temp
                    max_opt = (j, i)

            opt[j].append(temp)
    # Finish opt table ---------------------------

    result["score"] = max_score

    # TODO Need to figure out cigar string
    blah = []
    cigar = []
    prev = ""
    while max_opt[0] > 0:
        if d_opt[max_opt[0]][max_opt[1]] == 1:
            max_opt = (max_opt[0], max_opt[1] - 1)
            blah.insert(0, 'D')
        elif d_opt[max_opt[0]][max_opt[1]] == 0:
            max_opt = (max_opt[0] - 1, max_opt[1] - 1)
            blah.insert(0, 'M')
        elif d_opt[max_opt[0]][max_opt[1]] == -1:
            max_opt = (max_opt[0] - 1, max_opt[1])
            blah.insert(0, 'I')
        else:
            print("ruh Roh")

    cigar = []
    result["position"] = ref_pos - max_opt[1]
    result["cigar"] = "".join(blah)

    return result


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
    curr = 1
    count_table = {'$': 0}
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
# "index" "2019-nCoV.fa" "ref_index"
# "align" "ref_index" "reads_s1000.fa" "alignments"
