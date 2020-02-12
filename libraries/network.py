#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from sklearn import preprocessing
import networkx as nx
import matplotlib.pyplot as plt
import numpy 
import sys
from math import exp
import pickle as pk
import json
from os import path

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
    x = df['eval'].values
    listForces = []
    min = numpy.amin(x)
    max = numpy.amax(x)
    print(min,max)
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

def loadData(inputf):
    df = pd.read_csv(inputf, sep='\t', header=None)
    df.columns = ['node1','node2','eval']
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
    if len(df) == 0:
        print("NO DATA - QUITTING")
        exit()
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
    df_new = df.groupby(df['node1','node2']).aggregate(aggregation_functions)
    print("len of the dataframe after cleaning: " + str(len(df)))
    if len(df) == 0:
        print("NO DATA - QUITTING")
        exit()
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
        if node1 not in nodeList:
          nodeList.append(node1)
        if node2 not in nodeList:
          nodeList.append(node2)
        edgeList.append((node1,node2,eval))
    for n in nodeList:
          G.add_node(n)
    G.add_weighted_edges_from(edgeList)
    print("Loading NetworkX Data  - Finish")
    return G


def generateNetwork(inputf):
    expname = inputf.split(".")[0]
    df = loadData(inputf)
    df = removeSelfHit(df)
    df = removeDuplicates(df)
    #df = mergeHits(df)
    df = evalToForce(df)
    df = normalize(df)
    G = networkxLoad(df)
    # write the GML file for Gephi
    of = path.abspath(inputf).split(".")[0] 
    nx.write_gml(G,"{}.gml".format(of))
    # write json file for Cytoscape
    fjson = open("{}.cytoscape".format(of),"w")
    fjson.write(convert2cytoscapeJSON(G))
    fjson.close()
    
    
