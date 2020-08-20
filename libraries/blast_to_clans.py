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
from os import path

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

def generateClans(inputfile,use_pval):
    pvalue="1.0E-15"
    f = open(inputfile, 'rt')
    reader = csv.reader(f, delimiter='\t')
    listNamePos = []
    listName = []
    list_contacts = []
    count = 0
    textContacts = ""
    textSeq = ""

    print("creation of list of sequence names.")
    for row in reader:
        if row[0] not in listName:
            listName.append(row[0])
            listNamePos.append([row[0],count])
            textSeq += ">{}\nX\n".format(row[0])
            count += 1
    print("Done name list creation")

    print("creation  connections")
    f = open(inputfile, 'rt')
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
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

    of = path.abspath(inputfile).split(".")[0]
    fn = open("{}.clans".format(of),"w")
    fn.write(text.format(count,pvalue,textSeq,textContacts))
    fn.close()

    print("DONE - Generating CLANS file")
