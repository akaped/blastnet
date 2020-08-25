import argparse
from shutil import copyfile
from os import path, mkdir, system
from libraries.network import generateNetwork# import the required files to generate a cytoscape file that can be imported in gephi
from libraries.blast_to_clans import generateClans
from libraries.seqcounter import seqcounter
import subprocess

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

""")

def checkNCBIblastversion():
    version = str(subprocess.check_output(['blastn', '-version']))[10:16]
    print(f"You are using v{version} of NCBI Blast package")
    v1,v2,v3 = version.split(".")
    supportedv1 = ["2"]
    supportedv2 = ["8","9","10"]
    if v1 not in supportedv1 or v2 not in supportedv2:
        print(f"Your version of NCBI Blast {version}+ is not supported")
        print("The following versions are compatible with this script:")
        for i in supported:
            print(i)
        print("The script will execute anyway but it may generate execution errors.")

def makedb(fp,fn,runp,bt):
    print("GENERATING DB - START")
    if bt == "blastn":
        cmd = "makeblastdb -in {0} -out {1}/{2}_{3}_db -dbtype nucl ".format(fp,runp,fn,bt)
    elif bt == "psiblast" or bt == "blastr":
        cmd = "makeblastdb -in {0} -out {1}/{2}_{3}_db -dbtype prot ".format(fp,runp,fn,bt)
    elif bt == "parnassus":
        cmd = "makeblastdb -in {0} -out {1}/{2}_{3}_db -dbtype prot -blastdb_version 4".format(fp,runp,fn,bt)
    print(f"generating the database with this cmd: {cmd}")
    system(cmd)

    print("GENERATING DB - END\n")


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
    evalue = "0.1"
    cpu = "1"
    magicletter = ""
    blastype = ""
    matrix = "BLOSUM62"
    max_target_seqs = "1000"
    num_iterations = "1"
    inclusion_ethresh = "1"
    comp_based_stats = "2"

    if args.ifile:
        fp = path.abspath(args.ifile) #set full file path for input file
        fn = path.basename(fp).split(".")[0] # get filename
    else:
        print("No input file")
        exit()

    """ set blast arguments if specified otherwise default will be used """
    if args.evalue:
        evalue = args.evalue
    if args.cpu:
        cpu = args.cpu
    if args.matrix: #ok
        matrix = args.matrix
    if args.max_target_seqs:
        max_target_seqs = args.max_target_seqs
    if args.num_iterations:
        num_iterations = args.num_iterations
    if args.inclusion_ethresh and int(num_iterations) > 1:
        inclusion_ethresh =  args.inclusion_ethresh
    elif int(num_iterations) > 1:
        print("You need to specify the -inclusion_ethresh parameter if you want to run psiblast with more than one iteration")
        exit()
    if args.comp_based_stats:
        comp_based_stats = args.comp_based_stats
    else:
        if args.parnassus:
            comp_based_stats = "0"


    """ Check if the user inserted only the right arguments to run a version of blast, otherwise interrupt """
    if args.parnassus and args.n:
        print("You cant use the n flag while using parnassus")
        exit()
    if args.n and args.p:
        print("You can't run the script setting two flags, choose if you want to run blastp (-p) or blastn (-n)")
        exit()
    elif not args.n and not args.p and not args.parnassus and not args.counter:
        print("You should set the -p or -n flag, choose if you want to run blastp (-p) or blastn (-n)")
        exit()
    elif args.n and not args.p:
        print("* Selected search: NUCLEOTIDE SEARCH - evalue:{}".format(evalue))
        magicletter = "n"
        blastype = "blastn"
        # -S 1 Sets blast to search only forward sequences and not both ( as + reverse )
        cmd = 'blastn -db {}/{}db -query {}  -outfmt "6 qseqid qstart qend qlen qseq sseqid evalue pident bitscore sstart send slen length sseq" -out {}/{}_blastn.tsv -evalue {} -max_target_seqs ' + max_target_seqs + ' -num_threads {} -strand plus -max_hsps 1'
    elif not args.n and args.p and not args.parnassus:
        print("* Selected search: PROTEIN SEARCH - evalue:{}".format(evalue))
        magicletter = "p"
        blastype = "psiblast"
        cmd = 'psiblast -db {}/{}db -query {} -outfmt "6 qseqid qstart qend qlen qseq sseqid evalue pident bitscore sstart send slen length sseq" -out {}/{}_psiblast.tsv -evalue {} -matrix ' + matrix + ' -max_target_seqs ' + max_target_seqs + ' -num_iterations ' + num_iterations + ' -inclusion_ethresh ' + inclusion_ethresh + ' -num_threads {} -max_hsps 1 -comp_based_stats ' + comp_based_stats
    elif not args.n and args.parnassus:
        print("PARNASSUS binaries activated - using psiblast")
        print("* Selected search: PROTEIN SEARCH - evalue:{}".format(evalue))
        magicletter = "n"
        blastype = "parnassus"
        cmd = path.dirname(path.realpath(__file__)) + '/libraries/parnassus/psiblast -db {}/{}db -query {} -outfmt "6 qseqid qstart qend qlen qseq sseqid evalue pident bitscore sstart send slen length sseq" -out {}/{}_parnassus.tsv -evalue {} -matrix ' + matrix + ' -max_target_seqs ' + max_target_seqs + ' -num_iterations ' + num_iterations +' -inclusion_ethresh ' + inclusion_ethresh + ' -num_threads {} -max_hsps 1 -comp_based_stats ' + comp_based_stats
    elif not args.n and args.blastr:
        print("BLASTR binaries acrivated - using blastR")
        print("* Selected search: BLASTR search - evalue:{}".format(evalue))
        magicletter = "p"
        blastype = "blastr"
        cmd = path.dirname(path.realpath(__file__)) + '/libraries/blastR/scripts/blastallR.pl -p blastr -d {}/{}db -i {} -outfmt "6 qseqid qstart qend qlen qseq sseqid evalue pident bitscore sstart send slen length sseq" -o {}/{}_blastr.tsv -e {} -num_threads {} -max_hsps 1'
    if path.isfile(args.ifile):
        print("* Input file selected: {}".format(args.ifile))
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
        db = runp + "/{0}_{2}_db.{1}hr".format(fn,magicletter,blastype)#same here, shame on you #check if this is the file
        if path.isfile(db):
            print("* Database: DETECTED")
        elif args.n or args.p or args.parnassus:
            print("* Database: NOT DETECTED -- I will generate it")
            makedb(fp,fn,runp,blastype)
        else:
            print("I'm not generating the Db since it is not required")
        cmd = cmd.format(runp,fn,fp,rp,fn,evalue,cpu)
        print("* BLAST COMMAND - This is the blast command I'm running for you:")
        print(cmd)
        #print(cmd)
        return(fp,cmd,db,fn,runp,rp,evalue,blastype)
    else:
        print("Wrong path for your file, please check your input")
        exit()

def blastit(cmd,bt):
    try:
        print(f"* Running {bt.upper()}")
        system(cmd)
        print(f"* END {bt.upper()} RUN")
    except:
        print(f"! ERROR RUNNING {bt.upper()}")

def parsearg():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('ifile', action='store', type=str, help='your input file')
    parser.add_argument('-n', action="store_true",help="set the script to run BLASTN")
    parser.add_argument('-p', action="store_true",help="set the script to run BLASTP")
    parser.add_argument('-evalue', action="store", help="Cutoff Evalue of the blast search")
    parser.add_argument('-parnassus', action="store_true", help="uses alternate BLAST binaries for PARNASSUS analysis")
    parser.add_argument('-blastr', action="store_true", help="uses BlastR pearl scripts to run the analysis")
    parser.add_argument('-cpu', type=int, help='number of threads to run - good for cluster computers')
    parser.add_argument('-graphonly', action="store_true", help="The script will receive as input the tsv file and process only the graphs")
    parser.add_argument('-blastonly', action="store_true", help="The script will only generate the tsv file, but not process it")
    parser.add_argument('-counter', action="store_true", help="Counts number of seq per family and generates a csv file. necessary for bladerunner.py")
    parser.add_argument('-clans_use_eval', action="store_false", help="Normally pval is used to generate the CLANS output file, set this to switch to evalue")
    #parser.add_argument('-force_execution', action="store_true", help="Allows blastnet to be run with not supported/tested NCBI blast+")
    parser.add_argument('-matrix', action="store", help="Specify the type of matrix blast should use. Default: BLOSUM62")
    parser.add_argument('-max_target_seqs', action="store", help="Specify max number of sequences to retain per query. Default: 1000")
    parser.add_argument('-num_iterations', action="store", help="Number of iterations for psiblast. Default: 1")
    parser.add_argument('-inclusion_ethresh', action="store", help="Evalue of the inclusion threashold, this is required to use psiblast with more than one iteration")
    parser.add_argument('-comp_based_stats', action="store", help="Use composition-based statistics: Default 2")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    system("clear")
    banner()
    args = parsearg()
    #if not args.force_execution:
    checkNCBIblastversion()
    fp,cmd,db,fn,runp,rp,evalue,bt = checkinput(args)
    if args.counter:
        seqcounter(args.ifile,rp)
    if not args.graphonly:
        blastit(cmd,bt)
        bf = rp + "/" + fn + "_" + bt + ".tsv" #blast file
    if not path.isfile(rp + "/" + fn + "_" + bt + ".tsv"):
        print(f'The file {rp + "/" + fn + "_" + bt + ".tsv"} does\'t exist on this system')
        exit()

    if not args.blastonly and not args.counter:
        generateNetwork(bf)
        print("* GEPHI + CYTOSCAPE FILES GENERATED")
        if args.clans_use_eval:
            print("Using Eval to generate CLANS input file")
            generateClans(bf,args.ifile,False)
        else:
            print("Using pval to generate CLANS input file")
            generateClans(bf,args.ifile,True)
        print("* CLANS FILES GENERATED")
