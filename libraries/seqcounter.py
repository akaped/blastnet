from Bio import SeqIO

def seqcounter(ifile,rp):
    print("Counting seqs per family")
    families = []
    famcount = {}
    # extract from the seq_name the family name add it to a list
    for record in SeqIO.parse(ifile, "fasta"):
        fam =  str(record.id).strip().split("/")[0]
        if record.seq:
            families.append(fam)
    # count same elements in the list
    for fam in families:
        c = families.count(fam)
        # expected1 = int(pow(c,2))
        if fam not in famcount:
            famcount[fam]=c
    print("Finish counting seqs per family")
    if __name__ != "__main__":
        of = f"{rp}/seqfamcount.csv"
        csvtxt = ""
        for key,value in famcount.items():
            csvtxt += f"{key},{value}\n"
        fh = open(of, "w").write(csvtxt)
    #print(famcount)
    return famcount

if __name__ =="__main__":
  famcount = seqcounter(input("Please insert the path of the file you want to count"))
  print(famcount)
