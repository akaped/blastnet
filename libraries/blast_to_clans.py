#!/usr/bin/python3

""" Custom version of Clanstix, originally developed by M.Magnus """
""" Modified by P.Boccaletto """

""" Converts BLAST results to CLANS input format """
""" This script takes as input a csv file formatted with tabs, also known as tsv file, this file must have col1 col2 col3  """
""" Col1 = HIT, col2= query, col3 = e-value"""
""" Converts e-values to p-values """

""" conversion formula = Pvalue = 1-e^(-Evalue )"""
""" pay attention ! pvalue and evalue are almost identical when evalue < 0.01 !!! """

import csv
import sys
from math import exp
from os import path, listdir
from Bio import SeqIO
from sys import argv

#use_pval = True;
# set on false if you want to use eval instead
# maybe clans needs pval! So set it to true!

text = """
sequences={0}
<param>
maxmove=0.1
pval=1
usescval=false
complexatt=true
cooling=1.0
currcool=1.0
attfactor=10.0
attvalpow=1
repfactor=10.0
repvalpow=1
dampening=1.0
minattract=1.0
cluster2d=false
blastpath=''
formatdbpath=''
showinfo=false
zoom=1.0
dotsize=3
ovalsize=10
groupsize=4
usefoldchange=false
avgfoldchange=false
colorcutoffs=0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;
colorarr=(230;230;230):(207;207;207):(184;184;184):(161;161;161):(138;138;138):(115;115;115):(92;92;92):(69;69;69):(46;46;46):(23;23;23):
</param>
<seq>
{2}
</seq>

<hsp>
{3}
</hsp>
"""

def cal_pvalue(evalue):
    pvalue = 1-(1/exp(float(evalue)))
    return pvalue

def generateClans(blastdir,ifile,rp,fn,use_pval):
    pvalue="1.0E-15"
    #parsing the input fasta file into a dictionary <- this is necessary to get sequences
    seq_dict={}
    for record in SeqIO.parse(ifile,'fasta'):
        seq_dict[record.id]=str(record.seq)

    listNamePos = []
    listName = []
    list_contacts = []
    count = 0
    textContacts = ""
    textSeq = ""

    print("creation of list of sequence names.")
    for file in listdir(blastdir):
        f = open(path.join(blastdir,file), 'rt')
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            if row:
                if row[0] not in listName and not row[0] == "Search has CONVERGED!":
                    listName.append(row[0])
                    listNamePos.append([row[0],count])
                    textSeq += f">{row[0]}\n{seq_dict[row[0]]}\n"
                    count += 1
    #print(textSeq)
    print("Done name list creation")

    print("creation  connections")
    for file in listdir(blastdir):
        f = open(path.join(blastdir,file), 'rt')
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row:
                if not row[0] == "Search has CONVERGED!":
                    pos1 = 0
                    pos2 = 0
                    for r in listNamePos:
                        if r[0] == row[0]:
                            pos1 = r[1]
                        if r[0] == row[5]:
                            pos2 = r[1]
                    if use_pval: #decide if want to use pval or eval by setting var use_pval
                        textContacts += "{0} {1}:{2}\n".format(pos1,pos2,cal_pvalue(row[6]))
                    else:
                        textContacts += "{0} {1}:{2}\n".format(pos1,pos2,row[6])
        #print(textContacts)
    print("done creating connections")
    of = path.join(rp,fn)
    fn = open("{}.clans".format(of),"w")
    fn.write(text.format(count,pvalue,textSeq,textContacts))
    fn.close()

    print("DONE - Generating CLANS file")

if __name__ == "__main__":
    try:
        blastdir = argv[1] # dir where the results are stored
        fastafile = argv[2] # original input file
        rp = argv[3] # dir path for output folder
        fn = argv[4] # file name for output file
    except:
        print("The script must be used in this way")
        print("python3 blast_to_clans.py blastoutput.tsv queries.fasta")
        print("queries.fasta is the original file you used to generate the blast output")
        print("the output of the script is the same directory of your blastoutput.tsv file")
        exit()
    generateClans(blastdir,fastafile,rp,fn,True)
