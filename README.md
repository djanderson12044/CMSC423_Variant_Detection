# CMSC423_Variant_Detection
## Notes
The program assumes that all data is going to be in a "data" subdirectory of the current directory,
so I made it so that you don't have to put "data/" when referencing files in data. The 
program does it for you. (however, you still have to put the file extensions for files
other than the "ref_index". The "ref_index"'s file extension is added dynamically
in the code.)

Also, I was unable to finish the project and had a little bit of trouble
getting it to run via basic python command line functions (I had been
using the pycharm editor the entire time and its a bit different from
using like a straight up linux bash) (Everything will be in the writeup)

# Running The Program
This is just what worked for me when I was trying to run my program
on my WSL. Please let me know at
djanderson12044@gmail.com if you are having any trouble.
### Importing files
    python[version] -m pip install biopython
    python[version] -m pip install pysam 
These are the only files that I believe need to be imported (although
the pysam module isn't really used, it is still imported
and will error without it), but you
may need to try installing (in the same way as above):

    gzip
    pickle
    
### Running the program
Run

    chmod a+x fmmap.py
Then, to call the function

    python3 fmmap.py "function" "arguments".....
    
# Testing
### Gzipped reads
So the program currently gets its queries (reads) from a non-gziped fasta
file. However, if you wish to test with reads from a gziped file,
include line 120 in the code and indent accordingly.

### Output
Output still goes to the file specified in the align function call
but it is just a text file.

One issue is that each write is ended by "\r\n" which is the new line
for windows DOS stuff. So you might have to change it around a bit if you're on
a different os.
    