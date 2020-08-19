#!/usr/bin/env python
# coding: utf-8

#takes as input a csv or tsv file, and a evalue cutoff, loads the data in pandas and fiters the dataframe by this value.
#writes a clean csv/tsv file. If imported as library is able to take as input a pandas df and return a clean pandas df
import pandas as pd
from sys import argv

def loadData(file):
    if file.endswith(".csv"):
        df = pd.read_csv(file, sep=',', header=None)
    elif file.endswith(".tsv"):
        df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['node1','node2','eval']
    print(df)
    return df

def filterevalue(df,eval):
    print(f"Filtering by Evalue {str(eval)} - Start")
    to_drop = []
    for i in range(len(df)):
        node1 = df.iloc[i]['node1']
        node2 = df.iloc[i]['node2']
        evalu = df.iloc[i]['eval']
        #print(float(evalu),eval)
        if float(evalu) > eval:
          to_drop.append(i)
    df = df.drop(df.index[to_drop])
    print("len of the dataframe after cleaning: " + str(len(df)))
    return df

if __name__ == "__main__":
    if not argv[1] or not argv[0]:
        """
        This script takes as input a csv/tsv file and a evalue value.
        Please run it in this format: python3 filterbyevalue.py [file] [evalue]
        Example: python3 filterbyevalue.py myfile.csv 1e-10
        """
        exit()
    file = argv[1]
    file_name = file.split(".")[0]
    evalue = float(argv[2])
    df = loadData(file)
    result = filterevalue(df,evalue)
    df.to_csv(f'{file_name}_filtered_{str(evalue)}.csv', index = False)
