# Name: Tan Zi Ching
# ID: 1181102378
# Date: 25/9/2021
# Project

'''
In this project, you are required to develop a program. The program should:
1) Accept a gene sequence.
2) Perform BLAST to search for protein homologs.
3) Retrieve top 10 homologs.
4) Perform Multiple Sequence Alignment on to the sequences.
5) Extract the positon and amino acid of non-conserved regions.
'''

from Bio.Seq import *
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess # to run bash script

# blast function
def blast_function(query):
    query_file = open(query)
    query_record = SeqIO.read(query_file, format='fasta')
    results_handle = NCBIWWW.qblast(
                     'blastp',
                     'nr',
                     query_record.seq)
    res = results_handle.read()
    save_file = open("output/qblast_blastp.xml", 'w')
    save_file.write(res)
    save_file.close()
    results_handle = open("output/qblast_blastp.xml", 'r')
    blast_records = NCBIXML.parse(results_handle)

    E_VALUE_THRESH = 0.04
    count = 0
    hitlist = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            count = count + 1
            if count > 10: break
            for hsp in alignment.hsps:                
                hitrecord = SeqRecord(Seq(hsp.sbjct),
                                      id=str(count),
                                      name=alignment.hit_id,
                                      description=alignment.title[:100])
                hitlist.append(hitrecord)
    SeqIO.write(hitlist, "output/hitlist.fasta", "fasta")
    return hitlist

# multiple sequence alignment function
def msa_function(ifile):
    in_file = ifile
    out_file = "output/aligned.fasta"

    # get bash command for running clustelo
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True, force=True)
    #print(clustalomega_cline)

    # run bash command from python
    aligning = subprocess.Popen(str(clustalomega_cline),
                                stdout=subprocess.PIPE)
    output, error = aligning.communicate()   
    #print(output)

# find positions of nonconserved regions
def ncr_function(ifile):
    protfile = ifile
    homologues = []
    homolog_file = open(protfile)
    homolog_record = SeqIO.parse(homolog_file, 'fasta')
    for homolog in homolog_record:
        #print(homolog.seq)
        homologues.append(homolog)
    tofile = []
    nonconserved = []
    c=0
    for hom in homologues:
        s2 = str(hom.seq)
        if(c==0): 
            q = hom
            print("Query  :", s2)
            x="Q: "+s2
            tofile.append(x)
        else: 
            print("Subject:", s2)
            x="S: "+s2
            tofile.append(x)
            for i in range(len(q)):
                if(q[i] != s2[i]):
                    nonconserved.append(i)
        c+=1
    nonconserved.sort()
    #print(nonconserved)
    sym = []
    print("         ", end='')
    sym.append("   ")
    for i in range(len(q)):
        if i in nonconserved:
            print(",", end='')
            sym.append(",")
        else:
            print("*", end='')
            sym.append("*")

    sym = ''.join(sym)
    tofile.append(sym)

    with open("output/nonconserved-regions.txt", 'w') as f:
        f.writelines('\n'.join(tofile)+'\n')

# get DNA sequence

# user input
#seqinput = input("Please enter a DNA sequence: ")
 
# read from FASTA
valid = True
file = input("Please enter the FASTA file name (HBA2.fasta): ")
seqfile = open(file, 'r')
next(seqfile)
lines = []
lines = seqfile.readlines()
seqinput = ''.join(lines)
seqinput = seqinput.replace('\n', '')
seqinput = seqinput.upper()
seqinput = seqinput[0:len(seqinput)-len(seqinput)%3]
for base in seqinput:    # check for ambiguous bases
    if base not in ['A', 'T', 'G', 'C']:
        valid = False

dna = SeqIO.read(file, 'fasta')
inputid = dna.id
inputname = dna.name
inputdes = dna.description

if valid == False:
    print("Sequence not recognised!")
else:
    seq = Seq(seqinput)
    prot_seq = seq.translate()
    record = SeqRecord(prot_seq,
                        id = inputid,
                        name = inputname,
                        description = inputdes)
    filehandle = open('protquery.fasta', 'w')
    SeqIO.write(record, filehandle, 'fasta')
    filehandle.close()

    print("\nRunning BLAST...\n")
    hitlist = blast_function("protquery.fasta")
    print("Displaying top 10 homologues\n")
    for i in range(len(hitlist)):
        print(hitlist[i], end='\n')
    print("\n\nHit list generated as 'hitlist.fasta'\n")

    # create MSA file with query + hitlist
    iifile = []
    infilehandle = open('infile.fasta', 'w')

    iifile.append(record)

    for h in SeqIO.parse('output/hitlist.fasta', 'fasta'):
        iifile.append(h)

    SeqIO.write(iifile, 'infile.fasta', 'fasta')

    print("\nRunning Multiple Sequence Alignment...\n")
    msa_function("infile.fasta")
    print("\nMSA alignment generated as 'aligned.fasta'\n")

    print("\nDisplaying nonconserved regions\n")
    ncr_function("output/aligned.fasta")
    print("\n\nList of nonconserved regions generated as 'nonconserved-regions.fasta'\n")
    
    print("\n\n\n")