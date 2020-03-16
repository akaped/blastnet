import argparse
from shutil import copyfile
from os import path, mkdir, system
from libraries.network import generateNetwork # import the required files to generate a cytoscape file that can be imported in gephi
from libraries.blast_to_clans import generateClans

def banner():
    print("""

  ____   _              _    _   _        _
 |  _ \ | |            | |  | \ | |      | |
 | |_) || |  __ _  ___ | |_ |  \| |  ___ | |_
 |  _ < | | / _` |/ __|| __|| . ` | / _ \| __|
 | |_) || || (_| |\__ \| |_ | |\  ||  __/| |_
 |____/ |_| \__,_||___/ \__||_| \_| \___| \__|. . .

 An automated tool to perform Network Graph Analysis of proteic and nucleic sequences with BLAST.
 Author: Pietro Boccaletto 2020.
 Licence: MIT

 -----------------------------------------------------------------------------------------------

"""
"""
 BUGFIX - results, I advise in saving as blastp.tsv / blastn.tsv
 Otherwise overwritting of the files. TRY IT ON BLASTN



""")




def makedb(fp,fn,runp,magicletter):
    try:
        print("GENERATING DB - START")
        if magicletter == "n":
            cmd = "makeblastdb -in {0} -out {1}/{2}db -dbtype nucl ".format(fp,runp,fn)
        elif magicletter == "p":
            cmd = "makeblastdb -in {0} -out {1}/{2}db -dbtype prot ".format(fp,runp,fn)
        system(cmd)
        print("GENERATING DB - END\n")
    except:
        print("! ERROR WHILE MAKING DB")
        exit()

def checkinput(args):
    # runfiles
    # results
    fp = "" # File Path
    cmd = "" # Cmd to execute
    db = "" # DB file path
    fn = "" # Filename
    dp = "" # dir path
    rp = "" # result path
    runp = "" # runfiles path
    evalue = "10"
    cpu = 1

    if args.e:
        evalue = args.e
    if args.cpu:
        cpu = args.cpu
    if args.parnassus and args.n:
        print("You cant use the n flag while using parnassus")
        exit()
    if args.n and args.p:
        print("You can't run the script setting two flags, choose if you want to run blastp (-p) or blastn (-n)")
        exit()
    elif not args.n and not args.p and not args.parnassus:
        print("You should set the -p or -n flag, choose if you want to run blastp (-p) or blastn (-n)")
        exit()
    elif args.n and not args.p:
        print("* Selected search: NUCLEOTIDE SEARCH - evalue:{}".format(evalue))
        magicletter = "n"
        blastype = "blastn"
        cmd = 'blastn -db {}/{}db -query {}  -outfmt "6 qseqid sseqid evalue" -out {}/blastn.tsv -evalue {} -num_threads {}'
    elif not args.n and args.p and not args.parnassus:
        print("* Selected search: PROTEIN SEARCH - evalue:{}".format(evalue))
        magicletter = "p"
        blastype = "blastp"
        cmd = 'blastp -db {}/{}db -query {} -outfmt "6 qseqid sseqid evalue" -out {}/blastp.tsv -evalue {} -matrix BLOSUM62 -num_threads {}'
    elif not args.n and args.parnassus:
        print("PARNASSUS binaries activated - using blastp")
        print("* Selected search: PROTEIN SEARCH - evalue:{}".format(evalue))
        magicletter = "p"
        blastype = "blastparna"
        cmd = path.dirname(path.realpath(__file__)) + '/libraries/parnassus/blastp -db {}/{}db -query {} -outfmt "6 qseqid sseqid evalue" -out {}/blastparna.tsv -evalue {} -matrix BLOSUM62 -num_threads {}'
    if path.isfile(args.ifile):
        print("* Input file selected: {}".format(args.ifile))
        fp = path.abspath(args.ifile) #set full file path for input file
        fn = path.basename(fp).split(".")[0] # get filename
        if not path.isdir(fn): # if the directory doesnt exist creates it with name as filename
            mkdir(fn)
        dp = path.abspath(fn) # otherwise it exists and we can assign it to dp
        rp = dp + "/results"
        if not path.isdir(rp):
            mkdir(rp)
        runp = dp + "/runfiles"
        if not path.isdir(runp):
            mkdir(runp)
        copyfile(fp, runp + "/" + path.basename(fp))# this is not elegant and will not work in windows
        fp = runp + "/" + path.basename(fp) # this is not elegant and will not work in windows
        db = runp + "/{0}db.{1}hr".format(fn,magicletter)  #same here, shame on you #check if this is the file
        if path.isfile(db):
            print("* Database: DETECTED")
        else:
            print("* Database: NOT DETECTED -- I will generate it")
            makedb(fp,fn,runp,magicletter)
        cmd = cmd.format(runp,fn,fp,rp,evalue,cpu)
        fn = blastype
        #print(cmd)
        return(fp,cmd,db,fn,runp,rp,evalue)
    else:
        print("Wrong path for your file, please check your input")
        exit()

def blastit(cmd):
    try:
        print("* Running BLAST")
        system(cmd)
        print("* END Blast RUN")
    except:
        print("! ERROR RUNNING BLAST")

def parsearg():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('ifile', action='store', type=str, help='your input file')
    parser.add_argument('-n', action="store_true",help="set the script to run BLASTN")
    parser.add_argument('-p', action="store_true",help="set the script to run BLASTP")
    parser.add_argument('-e', action="store", type=float, help='evalue for blast run')
    parser.add_argument('-parnassus', action="store_true", help="uses alternate BLAST binaries for PARNASSUS analysis")
    parser.add_argument('-cpu', type=int, help='number of threads to run - good for cluster computers')
    parser.add_argument('-graphonly', action="store_true", help="The script will receive as input the tsv file and process only the graphs")
    parser.add_argument('-blastonly', action="store_true", help="The script will only generate the tsv file, but not process it")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    system("clear")
    banner()
    args = parsearg()
    if not args.graphonly:
        fp,cmd,db,fn,runp,rp,evalue = checkinput(args)
        blastit(cmd)
        bf = rp + "/" + fn + ".tsv" #blast file
    else:
        if path.isfile(args.ifile):
            bf = args.ifile
        else:
            print(f"The file {args.ifile} does't exist on this system")
            exit()
    if not args.blastonly:
        generateNetwork(bf)
        print("* GEPHI + CYTOSCAPE FILES GENERATED")
        generateClans(bf)
        print("* CLANS FILES GENERATED")
