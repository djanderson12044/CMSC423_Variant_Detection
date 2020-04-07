# Any and all necessary documentation is in the README (if there is any)

import sys
from Bio import SeqIO
import json


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
# ref_index: (.?)output file to store all the index information
#   thinking maybe Json
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

    return_index = {}  # What we are writing to ref_index using json

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

    # Adding to dictionary which is converted to json at the end of the program
    return_index["name"] = record.id
    return_index["length"] = len(text)
    return_index["sa"] = sa
    return_index["occ"] = occ
    return_index["lColumn"] = lColumn
    return_index["fColumn"] = fColumn

    with open(ref_index + ".json", "w") as outfile:
        json.dump(return_index, outfile)
        outfile.close()


# ref_index: (.?) File from the index definition
# reads: (.fa)
# alignments: (.sam)
def align(ref_index, reads, alignments):  # todo Finish implementing align method
    print("This is align Speakings")


# Helper Methods

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


# Calls the main() definition
main()
