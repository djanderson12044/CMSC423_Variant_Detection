# Any and all necessary documentation is in the README (if there is any)

import sys
from Bio import SeqIO

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
    record = SeqIO.read("data/" + ref, "fasta")
    print("ID: " + record.id)
    print()
    print(record.seq)
# ref_index: (.?) File from the index definition
# reads: (.fa)
# alignments: (.sam)
def align(ref_index, reads, alignments): # todo Finish implementing align method
    print("This is align Speakings")

# Calls the main() definition
main()