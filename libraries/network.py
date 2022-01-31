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
    print("Converting Evalues to Forces")
    x = df['eval'].values
    listForces = []
    min = numpy.amin(x)
    max = numpy.amax(x)
    print(f"Min Evalue {min}, Max Evalue {max}")
    for eval in x:
        force = max - (eval - min)
        listForces.append(force)
    df['eval'] = listForces
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
    print("Removing Self Hits - Start")
    """ this function removes the BLAST self hits """
    to_drop = []
    double_nodes = {}
    for i in range(len(df)):
        node1 = df.iloc[i]['node1']
        node2 = df.iloc[i]['node2']
        evalu = df.iloc[i]['eval']
        #print(node1,node2)
        if node1 == node2:
          to_drop.append(i)

    df = df.drop(df.index[to_drop])
    print("len of the dataframe after cleaning: " + str(len(df)))
    print("Removing Self Hits - Finish")
    return(df)

def removeDuplicates(df):
    print("Removing Duplicates - Start")
    df.sort_values('node1')
    df.drop_duplicates(keep=False,inplace=True)
    print("len of the dataframe after cleaning: " + str(len(df)))
    if len(df) == 0:
        print("NO DATA - QUITTING")
        exit()
    print("Removing Duplicates - Finish")
    return(df)

def mergeHits(df):
    print("Merging values for same node1 and node2 - start")
    print(df)
    aggregation_functions = {'eval': 'median',}
    df_new = df.groupby(df['node1', 'node2']).aggregate(aggregation_functions)
    print("len of the dataframe after merging: " + str(len(df)))
    print("Merging - Finish")
    return(df_new)

def normalize(df):
    print("Normalisation - Start")
    """ Normalisation of the column eval """
    x = df['eval']#.values.astype(float)
    x = x.to_frame()
    #scaler = preprocessing.StandardScaler()
    scaler = preprocessing.MinMaxScaler()
    x_scaled = scaler.fit_transform(x)
    df_normalized = pd.DataFrame(x_scaled)
    dflist = df_normalized.values
    df['eval'] = dflist
    print("Normalisation - Finish")
    return df


def networkxLoad(df):
    print("Loading NetworkX Data  - Start")
    G = nx.Graph()
    #get unique nodes, add edges and nodes in networkx
    nodeList = []
    edgeList = []
    for i in range(len(df)):
        node1 = df.iloc[i]['node1']
        node2 = df.iloc[i]['node2']
        eval = df.iloc[i]['eval']
        if eval == 0:
            eval = 0.001 # Avoid removal of edges that had the lowest evalue
        if node1 not in nodeList:
            nodeList.append(node1)
        if node2 not in nodeList:
            nodeList.append(node2)
        edgeList.append((node1, node2, eval))
    for n in nodeList:
        G.add_node(n)
    G.add_weighted_edges_from(edgeList)
    print("Loading NetworkX Data  - Finish")
    return G


def addSuperFamtoNetwork(G, file):
    with open(file, "r") as f:
        lines = f.readlines()
    fam_superfam = {}
    for line in lines:
        fam, superfam = line.split("\t")
        fam_superfam[fam] = superfam
    for node in G.nodes():
        # fallira per tutti quelli che non hanno un superfam
        node_fam = node.strip().split("__")[0]
        #print("NODE FAM", node_fam)
        try:
            G.nodes[node]['superfamily'] = fam_superfam[node_fam]
        except KeyError:
            #print(f'{node_fam} not in superfamlist')
            pass
    return G


def generateNetwork(input, evalue=0, superfam_file_path=""):
    #expname = input.split(".")[0]
    if path.isfile(input):
        df = loadData(input)
    elif path.isdir(input):
        df = loadDir(input)
    else:
        print("input not recognized")
        exit()
    if evalue != 0:
        df = filterevalue(df, evalue)
    df = removeSelfHit(df)
    if len(df) == 0:
        print("Lenght of the dataset after cleaning = 0 - I can't generate the gephi/cytoscape network")
        pass
    else:
        df = removeDuplicates(df)
    if len(df) == 0:
        print("Lenght of the dataset after cleaning = 0 - I can't generate the gephi/cytoscape network")
        pass
    else:
        #df = mergeHits(df) can't be turned on after removing the duplicates
        df = evalToForce(df)
        df = normalize(df)
        G = networkxLoad(df)
        if superfam_file_path:
            G = addSuperFamtoNetwork(G, superfam_file_path)
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
