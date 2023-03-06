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

#################
####Define constants and functions related to them    
#################

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
               'missing': 'black'
              }
def get_col(ref,alt):
    if ref==alt:
        return label_color["ref"]
    elif alt=="N":
        return label_color["missing"]
    elif alt=="-":
        return label_color["del"]
    else:
        return label_color[ref+">"+alt]


genes={}
genes[266]=('ORF1a',13460)
genes[13459]=('ORF1b',21555)
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

translationtable = {
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}
#source https://www.geeksforgeeks.org/dna-protein-python-3/ 
def translate(seq):
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= translationtable[codon]
    return protein

    
with open('libs/NC_045512.2.fasta', 'r') as file:
    refseq= file.read().partition("\n")[2]
    
    


def getSubstitutions(t):
    # Scan the table and get all mutations>50%
    pp=getPositionAndState(t,50,100,False)
    #build a list of annotation for each position
    finallist=[]
    for n,alt in pp:
        start=[i for i in genes if i<=n and genes[i][1]>=n]
        if start==[]:
            finallist+=[(str(n),"")]
        elif alt!="-":
            start=start[-1]
            #g=[genes[i][0] for i in genes if i<n and genes[i][1]>n][0]
            if start==266 and n-start>13460: #frameshift of ORF1b
                start-=1
            posCodonInProt=int((n-start)/3)+1
            posCodonInNuc=n-(n-start)%3
            posInCodon=(n-start)%3
            refcodon=refseq[posCodonInNuc-1:posCodonInNuc+2]
            newcodon = refcodon[:posInCodon] + alt + refcodon[posInCodon + 1:]
            for n2 in range(n-posInCodon,n-posInCodon+3):
                if n2!=n and n2 in [i[0] for i in pp]:
                    alt2=[i[1] for i in pp if i[0]==n2][0]
                    posInCodon2=posInCodon+n2-n
                    newcodon = newcodon[:posInCodon2] + alt2 + newcodon[posInCodon2 + 1:]
            AAinit=translate(refcodon)
            AAnew=translate(newcodon)
            if AAinit!=AAnew:
                if posCodonInProt==1:
                    AAnew="?"
                finallist+=[(str(n),AAinit+str(posCodonInProt)+AAnew)]
    return(finallist)
            
       
            
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
    columns=list(set(t.columns)-set(['POS', 'REF']))
    return(t[columns].sum(axis=1))

# Loop through a table and keep only the positions where an alternative allele represent
# more than min% of the total number of samples
def getPositionAndState(t,percentmin,percentmax,addmissing):
    statelist=["A","C","G","T","-"]
    if addmissing : statelist+=["N"]
    #absolute minimum number of sample
    columns=list(set(t.columns)-set(['POS', 'REF']))
    n_min=sumperline(t)/100*percentmin
    n_max=sumperline(t)/100*percentmax
    # for each possible reference allele :
    totalposlist=[]
    for s in statelist:
        p=t[(t["REF"]!=s) & (s!=0) & (t[s]>=n_min) & (t[s]<=n_max)]["POS"]
        totalposlist+=[(i,s) for i in p]
    return totalposlist


# Loop through the tables and keep only the positions wher an alternative allele represent
# more than min% of the total number of samples
def getpositions(tablelist,percentmin=0,percentmax=100,addmissing=False):
    totalposlist=[]
    for t in tablelist:
        totalposlist+=getPositionAndState(t,percentmin,percentmax,addmissing)
    totalposlist=list(set([i[0] for i in totalposlist]))
    totalposlist.sort()
    return totalposlist

    

def autolabel(ax,x,text,textcolor):
    ax.text(x, 50, text, color=textcolor,ha='center', va='center', rotation=90,weight="bold")

    
def sumperline(t):
    return(t[list(set(t.columns)-set(['POS', 'REF']))].sum(axis=1))
    
def addpos(ref,alt,all_bottoms,axx,all_pos_toplot,x_names):
    values=all_pos_toplot[alt].copy()
    if ref in ["A","C","G","T"]:
        values[all_pos_toplot['REF'] != ref] = 0
    values=values/sumperline(all_pos_toplot)*100
    axx.bar(x_names, values,
            bottom=all_bottoms,
            color=get_col(ref,alt),
            edgecolor="none",width=1)
    return values
    
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


def bighist(tablelist,poslist,y_names,mytitle="",suptables=[],PDFname=""):
    
    fig = plt.figure(figsize=(len(poslist)/3+5,2*(len(tablelist)+len(suptables))), constrained_layout=False)
    ax = fig.add_gridspec(nrows=len(tablelist)+len(suptables)+1, ncols=1, hspace=0).subplots(sharex=True)
    x_names=[str(i) for i in poslist]
    for i in range(len(suptables)):
        ax[i].bar(x_names , suptables[i] , color="grey" , edgecolor="none" , width=1)
    for i,t in enumerate(tablelist):
        axx=ax[i+len(suptables)]
        axx.set(yticks=[25,50,75])
        axx.set_ylim((0,100))
        axx.margins(0, 0)  
        nb_sample=sum(t.iloc[0][list(set(t.columns)-set(['POS', 'REF']))])
        if nb_sample<10000:str_nb_sample="\nn="+str(nb_sample)+""
        else: str_nb_sample="\nn="+str(int(nb_sample/1000))+"K"
        axx.set_ylabel(y_names[i]+str_nb_sample, fontsize=15)
        
        all_pos_toplot=t.loc[t["POS"].isin(poslist)]
        all_bottoms=all_pos_toplot['A']-all_pos_toplot['A']
        all_bottoms+=addpos("*","-",all_bottoms,axx,all_pos_toplot,x_names)
        #loop in all the 12 combinaisons:
        for ref in ["A","C","G","T"]:
            for alt in ["A","C","G","T"]:
                if alt!=ref:
                    all_bottoms+=addpos(ref,alt,all_bottoms,axx,all_pos_toplot,x_names)
            all_bottoms+=addpos(ref,ref,all_bottoms,axx,all_pos_toplot,x_names)
            addpos(ref,"N",all_bottoms,axx,all_pos_toplot,x_names)
        for p,lab in getSubstitutions(t):
            if p in x_names: autolabel(axx,p,lab,"black")
    ax_last=ax[len(tablelist)+len(suptables)]
    addgenenames(ax_last,x_names)
    ax_last.tick_params(axis='x',bottom=True,labelrotation=90, labelsize=20)
    legend=[mpatches.Patch(color=label_color[i], label=i.replace("T","U")) for i in label_color]
    fig.legend(handles=legend,loc='upper right',title=mytitle)
    fig.subplots_adjust(right=0.9)
    if PDFname!="":
        fig.savefig(PDFname, bbox_inches='tight')
        

        

        
def createsuptable(positions,filename, colname):
    t=pd.read_csv(filename,sep="\t")
    return t[t['pos_in_genome'].isin(poslist)][colname]


