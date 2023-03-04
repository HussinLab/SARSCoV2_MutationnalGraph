#!/usr/bin/env python3
import os
from os import path
import sys
import csv

import numpy as np
import pandas as pd


from functools import reduce
import math

#graphic libraries
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.cm import hsv
from matplotlib.lines import Line2D
import seaborn as sns
sns.set_theme(style="ticks", color_codes=True)


def autolabel(ax,x,text,textcolor):
    ax.text(x, 50, text, color=textcolor,ha='center', va='center', rotation=90,weight="bold")
    
#dictionnary to label each mutation
nucsub_AAname = {}
def load_mut_names(file):
    for i,(nuc,AA) in pd.read_csv(file,sep="\t",header=None).iterrows():
        if AA[0] == AA[-1]: AA="" # do not name synonymous mutations
        nucsub_AAname[nuc]=AA


min_val_AAlabel=10 # % of alt alleles to add amino acide label
def def_min_val_label(n):
    global min_val_AAlabel
    min_val_AAlabel=n
        
#colors of the different mutations
label_color = {'A>C': 'darkorange',
               'A>G': 'gold',
               'A>T': 'mediumpurple',
               'C>A': 'limegreen',
               'C>G': 'violet' ,
               'C>T': 'dodgerblue',
               'G>A': 'mediumblue',
               'G>C': 'fuchsia',
               'G>T': 'forestgreen',
               'T>A': 'indigo',
               'T>C': 'dimgrey',
               'T>G': 'red',
               'ref': 'white',
               'del': 'silver',
               'missing': 'black'}


def get_col(ref,alt):
    if ref==alt:
        return label_color["ref"]
    elif alt=="N":
        return label_color["missing"]
    elif alt=="-":
        return label_color["del"]
    else:
        return label_color[ref+">"+alt]

    
with open('libs/NC_045512.2.fasta', 'r') as file:
    refseq= file.read().partition("\n")[2]
    
#source https://www.geeksforgeeks.org/dna-protein-python-3/ 
def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

with open('libs/NC_045512.2.fasta', 'r') as file:
    refseq= file.read().partition("\n")[2]

    
def getMut(n,alt):
    start=[i for i in genes if i<n and genes[i][1]>n]
    if start==[]:
        return("IG")
    else:
        start=start[0]
    #g=[genes[i][0] for i in genes if i<n and genes[i][1]>n][0]
    if start==266 and n-start>13460: #frameshift of ORF1b
        start-=1
    codonstart=n-(n-start)%3
    posincodon=(n-start)%3
    refcodon=refseq[codonstart-1:codonstart+2]
    newcodon = refcodon[:posincodon] + alt + refcodon[posincodon + 1:]
    if(translate(refcodon)!=translate(newcodon)):
        return(translate(refcodon)+str(int((n-start)/3)+1)+translate(newcodon))
    else:
        return("")

    
genes={}
genes[266]=('ORF1ab',21555)
genes[21563]=('Spike',25384)
genes[25393]=('ORF3a',26220)
genes[26245]=('E',26472)
genes[26523]=('M',27191)
genes[27202]=('ORF6',27387)
genes[27394]=('ORF7a',27759)
genes[27756]=('ORF7b',27887)
genes[27894]=('ORF8',28259)
genes[28274]=('N',29533)
genes[29558]=('ORF10',29674)

            
def addgenenames(ax, x_names):
    n_pos=len(x_names)
    ax.set(yticks=[25,50,75])
    ax.set_ylabel("gene\nnames",  fontsize=20)
    ax.set_ylim((0,100))
    ax.margins(0, 0)  
    colcol="grey"
    for lower_bound in genes.keys():
        (gene_name,upper_bound)=genes[lower_bound]
        pospos=[i for i in x_names if (int(i)>lower_bound and int(i)<upper_bound)]
        if pospos!=[]:
            span=[(i in pospos)*100 for i in x_names]
            ax.bar(x_names , span , bottom=[0]*n_pos , color=colcol , edgecolor="none" , width=1)
            if colcol=="grey": colcol="black"
            else             :  colcol="grey"
            mean=0
            for i in range(n_pos):
                if list(x_names)[i] in pospos:
                    mean+=i/len(pospos)
            mean=int(mean)
            autolabel(ax,mean,gene_name,"white")
            
            
#Create a list of panda table containing all the numbers for each of the 29903 positions
def openfiles(inputfiles):
    tablelist=[]
    for file in inputfiles:
        if path.isfile(file):
            t=pd.read_csv(file,sep="\t")
            tablelist.append(t)
        else:
            print("ERROR NO SUCH FILE : "+file)
            return None;
    return tablelist

def sumperline(t):
    return(t[set(t.columns)-set(['POS', 'REF'])].sum(axis=1))

# Loop through the tables and keep only the positions wher an alternative allele represent
# more than min% of the total number of samples
def getpositions(tablelist,percentmin=0,percentmax=100,addmissing=False):
    totalposlist=[]
    statelist=["A","C","G","T","-"]
    if addmissing : statelist+=["N"]
    for t in tablelist:
        #absolute minimum number of sample
        n_min=t[set(t.columns)-set(['POS', 'REF'])].sum(axis=1)/100*percentmin
        n_max=t[set(t.columns)-set(['POS', 'REF'])].sum(axis=1)/100*percentmax
        # for each possible reference allele :
        for s in statelist:
            p=t[(t["REF"]!=s) & (s!=0) & (t[s]>=n_min) & (t[s]<=n_max)]["POS"]
            totalposlist+=list(p)#[str(i)+"."+s for i in p]
    totalposlist=list(set(totalposlist))
    totalposlist.sort()
    return totalposlist
    




    
def bighist(tablelist,poslist,y_names,mytitle="",suptables=[],PDFname=""):
    
    fig = plt.figure(figsize=(len(poslist)/3+5,2*(len(tablelist)+len(suptables))), constrained_layout=False)
    ax = fig.add_gridspec(nrows=len(tablelist)+len(suptables)+1, ncols=1, hspace=0).subplots(sharex=True)
    x_names=[str(i) for i in poslist]
    n_table=len(tablelist)
    for i in range(len(suptables)):
        ax[i].bar(x_names , suptables[i] , color="grey" , edgecolor="none" , width=1)
    for i in range(n_table):
        axx=ax[i+len(suptables)]
        nb_sample=sum(tablelist[i].iloc[0][["A","C","G","T","-","N"]])
        if nb_sample>2:
            if nb_sample<10000:
                str_nb_sample="\nn="+str(nb_sample)+""
            else:
                str_nb_sample="\nn="+str(int(nb_sample/1000))+"K"
        else :
            str_nb_sample=""
        axx.set(yticks=[25,50,75])
        axx.set_ylabel(y_names[i]+str_nb_sample,  fontsize=15)
        axx.set_ylim((0,100))
        axx.margins(0, 0)  
        all_pos_toplot=tablelist[i].loc[tablelist[i]["POS"].isin(poslist)]
        all_bottoms=all_pos_toplot['A']-all_pos_toplot['A']
        def addpos(ref,alt,all_bottoms):
            values=all_pos_toplot[alt].copy()
            if ref in ["A","C","G","T"]:
                values[all_pos_toplot['REF'] != ref] = 0
            values=values/nb_sample*100
            axx.bar(x_names, values,
                    bottom=all_bottoms,
                    color=get_col(ref,alt),
                    edgecolor="none",width=1)
            for j in range(len(values)):
                if values.iloc[j]>min_val_AAlabel:
                    idx=str(all_pos_toplot.iloc[j]["POS"])+ref+">"+alt
                    if idx in nucsub_AAname:
                        if ref+">"+alt in ["T>A","G>A","T>C"]:
                            autolabel(axx,j,getMut(all_pos_toplot.iloc[j]["POS"],alt),"lightgrey")
                        else: autolabel(axx,j,getMut(all_pos_toplot.iloc[j]["POS"],alt),"black")
            return values
        all_bottoms+=addpos("*","-",all_bottoms)
        #loop in all the 12 combinaisons:
        for ref in ["A","C","G","T"]:
            for alt in ["A","C","G","T"]:
                if alt!=ref:
                    values=addpos(ref,alt,all_bottoms)
                    all_bottoms+=values
            values=addpos(ref,ref,all_bottoms)
            all_bottoms+=values
            values=addpos(ref,"N",all_bottoms)
    ax_last=ax[n_table+len(suptables)]
    addgenenames(ax_last,x_names)
    ax_last.tick_params(axis='x',bottom=True,labelrotation=90, labelsize=20)
    legend=[mpatches.Patch(color=label_color[i], label=i.replace("T","U")) for i in label_color]
    fig.legend(handles=legend,loc='center right',title=mytitle, bbox_to_anchor=(1.05, 0.5))
    fig.subplots_adjust(right=0.9)
    if PDFname!="":
        fig.savefig(PDFname, bbox_inches='tight')
        
def createsuptable(positions,filename, colname):
    t=pd.read_csv(filename,sep="\t")
    return t[t['pos_in_genome'].isin(poslist)][colname]


