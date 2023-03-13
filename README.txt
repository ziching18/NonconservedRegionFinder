Please run "NonconservedRegionFinder.py".

====================================================================

About the Program
-----------------

This program accepts a FASTA file from the user and does 3 things:
- run BLAST on the sequence,
- perform multiple sequence alignment on the top 10 homologues,
- extract the positions of non-conserved regions.

A FASTA file is provided for the users (HBA2.fasta) and has been written in the input message to prompt the user to enter the correct file name. If there are ambiguous codes in the sequence, like letters other than A, T, G and C, the program will terminate. Else, all 3 functions will run, and the required files will be generated in the 'output' folder. Contained in the folder "sample output" are what the output should be.

====================================================================

Prerequisite
------------

Please download the clustal omega executable from http://www.clustal.org/omega/
The main 'clustalo.exe' and the '.dll' files must be included within the directory.

====================================================================

About the Output
----------------

The program output are:
1. 'qblast_blastp.xml' - the top 50 homologues of the user input sequence generated from BLAST
2. 'hitlist.fasta' - the top 10 homologues extracted from the xml file
3. 'aligned.fasta' - the output of multiple sequence alignment
4. 'nonconserved-regions.txt' - the positions of nonconserved regions will be symbolised with a "," instead of a "*" underneath the sequences, the first sequence is the query sequence, followed by the rest of the subject sequences.

====================================================================

Other files
-----------

There are several files included in the folder that are meant for the multiple sequence alignment program. These include "clustalo.exe" and its dependent .dll files.

====================================================================

Thank you!
- Zi Ching