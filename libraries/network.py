#!/usr/bin/env python
# coding: utf-8

import ray
ray.init()
from os import path, listdir
import modin.pandas as pd
#import pandas as pd
from sklearn import preprocessing
import networkx as nx
import matplotlib.pyplot as plt
import numpy
import sys
from math import exp
import pickle as pk
import json
from sys import argv
from filterbyevalue import filterevalue

def convert2cytoscapeJSON(G):
    # load all nodes into nodes array
    final = {}
    final["nodes"] = []
    final["edges"] = []
    for node in G.nodes():
        nx = {}
        nx["data"] = {}
        nx["data"]["id"] = node
        nx["data"]["label"] = node
        final["nodes"].append(nx.copy())
    #load all edges to edges array
    for edge in G.edges():
        nx = {}
        nx["data"]={}
        nx["data"]["id"]=edge[0]+edge[1]
        nx["data"]["source"]=edge[0]
        nx["data"]["target"]=edge[1]
        final["edges"].append(nx)
    return json.dumps(final)

def evalToForce(df):
    print("## EVALUE TO FORCE ##")
    print("Converting Evalues to Forces")
    x = df['eval'].values
    listForces = []
    min = numpy.amin(x)
    max = numpy.amax(x)
    print(f"Min Force {min}, Max Force {max}")
    for eval in x:
        force = max - (eval - min)
        listForces.append(force)
    df['eval'] = listForces
    print('------------------------------')
    return df

# Obsolete, pval is not a good measure
"""def cal_pvalue_force(df):
    x = df['eval'].values
    listpval = []
    for i in x:
        pval = 1-(1/exp(i))
        force = pval - 1
        listpval.append(force)
    df['eval'] = listpval
    return df
"""

""" deprecated """
def loadData(inputf):
    if inputf.endswith(".csv"):
        df = pd.read_csv(inputf, sep=',', header=None)
    elif inputf.endswith(".tsv"):
        df = pd.read_csv(inputf, sep='\t', header=None)
    else:
        print(f"File {inputf} not valid")
        exit()
    if len(df.columns) == 14:
        df.columns = ["node1","qstart","qend","qlen","qseq","node2","eval","pident","bitscore","sstart","send","slen","length","sseq"]
        df = df.drop(columns=["qend","qstart","qlen","qseq","pident","bitscore","sstart","send","slen","length","sseq"])
    else:
        df.columns = ["node1","node2","eval"]
    print("CSV loaded into Pandas Dataframe - lenght of dataset : {}".format(str(len(df))))
    return df

def loadDir(input):
    df = pd.DataFrame()
    #print(input)
    for file in listdir(input):
        #print(input,file)
        if file.endswith(".csv"):
            tempdf = pd.read_csv(path.join(input,file), sep=',', header=None)
        elif file.endswith(".tsv"):
            tempdf = pd.read_csv(path.join(input,file), sep='\t', header=None)
        else:
            print(f"File {input} not valid")
            exit()
        df = df.append(tempdf)
    df.columns = ["node1","qstart","qend","qlen","qseq","node2","eval","pident","bitscore","sstart","send","slen","length","sseq"]
    df = df.drop(columns=["qend","qstart","qlen","qseq","pident","bitscore","sstart","send","slen","length","sseq"])
    print("CSV loaded into Pandas Dataframe - lenght of dataset : {}".format(str(len(df))))
    return df


def removeSelfHit(df):
    print("## REMOVE SELF HITS ##")
    """ this function removes the BLAST self hits - OK TESTED """
    to_drop = []
    for i in range(len(df)):
        node1 = df.iloc[i]['node1']
        node2 = df.iloc[i]['node2']
        if str(node1) == str(node2):
            to_drop.append(i)
    df = df.drop(df.index[to_drop])
    print("Lenght of the dataframe after REMOVE SELF HITS: " + str(len(df)))
    if len(df) == 0:
        print("ERROR: Lenght of the dataframe = 0 - I can't generate the gephi/cytoscape network")
        exit()
    print('------------------------------')
    return(df)


def removeDuplicates(df):
    print("## REMOVE DUPLICATES HITS ##")
    ''' Drop duplicates removes only! full line duplicates ...
    CHECKED - OK '''
    df.sort_values('node1')
    df.drop_duplicates(keep="first", inplace=True)
    print("Lenght of the dataframe after REMOVE DUPLICATE HITS: " + str(len(df)))
    if len(df) == 0:
        print("ERROR: Lenght of the dataframe = 0 - I can't generate the gephi/cytoscape network")
        exit()
    print('------------------------------')
    return(df)


def mergeHits(df):
    ''' CHECKED - OK '''
    print("## MERGE HITS ##")
    df = df.groupby(['node1','node2'], as_index=False)['eval'].mean()
    print("Lenght of the dataframe after MERGE HITS: " + str(len(df)))
    if len(df) == 0:
        print("ERROR: Lenght of the dataframe = 0 - I can't generate the gephi/cytoscape network")
        exit()
    print('------------------------------')
    return(df)


def normalize(df):
    ''' CHECKED - OK '''
    print("## NORMALIZING EDGES ##")
    """ Normalisation of the column eval """
    x = df['eval']#.values.astype(float)
    x = x.to_frame()
    #scaler = preprocessing.StandardScaler()
    scaler = preprocessing.MinMaxScaler()
    x_scaled = scaler.fit_transform(x)
    df_normalized = pd.DataFrame(x_scaled)
    dflist = df_normalized.values
    df['eval'] = dflist
    print('------------------------------')
    return df


def networkxLoad(df):
    ''' CHECKED - OK '''
    print("## Loading NetworkX DATA ##")
    G = nx.Graph()
    #get unique nodes, add edges and nodes in networkx
    nodeList = []
    edgeList = []
    min_force = numpy.amin(df['eval'].values)
    for i in range(len(df)):
        node1 = df.iloc[i]['node1']
        node2 = df.iloc[i]['node2']
        eval = df.iloc[i]['eval']
        if node1 not in nodeList:
            nodeList.append(node1)
        if node2 not in nodeList:
            nodeList.append(node2)
        if eval != min_force:
            ''' remove all the edges that have eval = min evalue '''
            edgeList.append((node1, node2, eval))
    for n in nodeList:
        G.add_node(n)
    G.add_weighted_edges_from(edgeList)
    print('------------------------------')
    return G


def addFamAndSuperFamtoNetwork(G, file):
    ''' CHECKED - OK '''
    ''' Allows to add the attribute Family and Superfamily to each node '''
    print("## ANNOTATING NODES ##")
    with open(file, "r") as f:
        lines = f.readlines()
    fam_superfam = {}
    for line in lines:
        line = line.split("\t")
        fam = line[0].strip()
        superfam = line[1].strip()
        fam_superfam[fam] = superfam
    for node in G.nodes():
        # fallira per tutti quelli che non hanno un superfam
        node_fam = node.strip().split("__")[0]
        #print("NODE FAM", node_fam)
        try:
            G.nodes[node]['superfamily'] = fam_superfam[node_fam]
        except KeyError:
            G.nodes[node]['superfamily'] = "not_assigned"
            pass
        try:
            G.nodes[node]['family'] = node_fam
        except Exception as e:
            print(e)
            pass
    print('------------------------------')
    return G


def removeBigEval(df):
    ''' CHECKED - OK '''
    print("## KEEP ONLY HIT WITH SMALLER EVAL ##")
    df.sort_values('eval')
    df.drop_duplicates(subset = ['node1', 'node2'], keep="first", inplace=True)
    print("Lenght of the dataframe after KEEP ONLY HIT WITH SMALLER EVAL: " + str(len(df)))
    if len(df) == 0:
        print("ERROR: Lenght of the dataframe = 0 - I can't generate the gephi/cytoscape network")
        exit()
    print('------------------------------')
    return(df)


def generateNetwork(input, evalue=0, superfam_file_path=""):
    #expname = input.split(".")[0]
    ''' Load of files and directories '''
    if path.isfile(input):
        df = loadData(input)
    elif path.isdir(input):
        df = loadDir(input)
    else:
        print("input not recognized")
        exit()
    ''' filtering by evalue '''
    if evalue != 0:
        df = filterevalue(df, evalue)
    ''' remove hits where node1 == node2 '''
    df = removeSelfHit(df)
    ''' remove hits that are duplicates, keeps only the first one '''
    df = removeDuplicates(df)
    #''' merge hits by evalue (mean) that have the same node1 node2 name '''
    #df = mergeHits(df)
    ''' keeps the hit with the lowest evalue discards the others '''
    df = removeBigEval(df)
    ''' Transforms to a force, where lower evalue (low value) corresponds
    to biggest attraction (big value)'''
    df = evalToForce(df)
    ''' normalizes the forces in the space from 0 to 1 '''
    df = normalize(df)
    ''' Generates the network, removed the edges that have force value = min normalized force value'''
    G = networkxLoad(df)
    if superfam_file_path:
        G = addFamAndSuperFamtoNetwork(G, superfam_file_path)
    print("FINAL STAGE --- Writing gml file. Please wait")
    # write the GML file for Gephi
    of = path.abspath(input).split(".")[0]
    nx.write_gml(G,f"{of}.gml")
    # write json file for Cytoscape
    #fjson = open("{}.cytoscape".format(of),"w")
    #fjson.write(convert2cytoscapeJSON(G))
    #fjson.close()


if __name__ == "__main__":
    try:
        superfam_file_path = argv[3]
    except Exception:
        superfam_file_path = ""
    generateNetwork(argv[1], float(argv[2]), superfam_file_path)
