#!/usr/bin/env python
# coding: utf-8

# In[2]:


import baltic as bt
import requests
from io import StringIO as sio
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon, PathPatch
from matplotlib.path import Path
from matplotlib.collections import LineCollection
from matplotlib import cm
import matplotlib.patheffects as path_effects
import numpy as np
from datetime import datetime as dt
import os
import PyAstronomy
from PyAstronomy import pyasl
import datetime
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon


from scipy import stats
import pandas as pd

def get_all_terminals(node):
    """
    Get all terminal nodes that descend from a given node
    """
    terminals = []
    for i in node.children:
        if i.branchType == "leaf":
            terminals.append(i)
        if i.branchType != "leaf":
            terminals.extend(get_all_terminals(i))
    return terminals


df=pd.read_csv("tmrcas.txt", sep="\t")
df2=pd.read_csv("Mpox_2poch_combined.log", sep="\t")
df['tmrca(ingroup)']=df2['tmrca(ingroup)'].to_list()[0:len(df)]

mrcas={}

for c in list(df.columns[1:]):
    vals=df[c].to_list()
    burin=int(round(len(vals)*0.1,0))
    vals=vals[burin:]
    mrcas[c]=[2023.4137-float(xx) for xx in vals]

mrcas['age(apobec3.transition)']=df2['age(apobec3.transition)'].to_list()[0:len(df)]


# In[33]:


treeFile="hmpxv.nex" 
ll=bt.loadNexus(treeFile) 
ll.traverse_tree()
ll.treeStats()



for i in ll.Objects:
    if i.branchType == "leaf":
        if "VSP" in i.name:
            i.traits["rate_signDistribution"]="NEW"

        elif "TRM" in i.name:
            i.traits["rate_signDistribution"]="NEW"
        else:
            i.traits["rate_signDistribution"]="NO"



fig, ax1 = plt.subplots(figsize=(8, 8))
x_attr=lambda k: k.absoluteTime 
c_func=lambda k: 'black' 
cmap=mpl.cm.viridis

                
x_attr3=[]
for i in ll.Objects:
    if i.branchType == "node":
        nds=[[2023.4137-x for x in i.traits['height_95%_HPD']], i.y]
        x_attr3.append([nd for nd in nds])
        
    
effects=[path_effects.Stroke(linewidth=2, foreground='white'),
                 path_effects.Stroke(linewidth=0.5, foreground='black')]
kwargs={'ha':'left','va':'center','size':8,'path_effects':effects}


text_func=lambda k: k.name
t_target=lambda k: k.branchType=='leaf'
n_target=lambda k: k.branchType=='node' and len(k.children)==1



def  c_func2(k):
    if "Cameroon" in k.name:
        cc='#FFA000'
    elif k.name in ["unpub|VSP191|Nigeria|Abia|2022-12-30", "unpub|VSP189|Nigeria|Abia|2023-01-01"]:

        cc="#00356B"
    elif  k.name in ["unpub|TRM288|Nigeria|Akwa-Ibom|2022-09-28", "unpub|VSP199|Nigeria|Akwa-Ibom|2022-01-01"]:
        cc="#4997D0"
    elif "1971" in k.name:
        cc="#00356B"
    else:
        cc="#89CFF0"
    return(cc)
        
for xx in range(2010, 2025, 1):
    plt.axvspan(xx, xx+1, facecolor='#3AA8C1', alpha=0.1, zorder=1, edgecolor=None)  


ll.plotTree(ax1,x_attr=x_attr,colour=c_func,zorder=100, width=1.2) 

ll.plotPoints(ax1,x_attr=x_attr,size=35,colour=c_func2,zorder=100) 


ax2_2 = ax1.twinx()
sns.kdeplot(mrcas['tMRCA(Abia_human)'], ax=ax2_2, color='#0070BB', shade=True, alpha=0.2,zorder=1)
ax1.scatter(np.median(mrcas['tMRCA(Abia_human)']),174,s=80,color="#0070BB",zorder=500, edgecolor='black')
ax2_2.set_ylim(0, 5)
ax1.plot([np.median(mrcas['tMRCA(Abia_human)']), np.median(mrcas['tMRCA(Abia_human)'])], [0, 175], color='#0070BB', linestyle='--', linewidth=0.5)


ax2_5 = ax1.twinx()
sns.kdeplot(mrcas['tmrca(ingroup)'], ax=ax2_5, color='#003262', shade=True, alpha=0.2,zorder=1)
ax1.scatter(np.median(mrcas['tmrca(ingroup)']),167,s=80,color="#003262",zorder=500, edgecolor='black')
ax2_5.set_ylim(0, 7)
ax1.plot([np.median(mrcas['tmrca(ingroup)']), np.median(mrcas['tmrca(ingroup)'])], [0, 167], color='#003262', linestyle='--', linewidth=0.5)


ax2_1 = ax1.twinx()
sns.kdeplot(mrcas['age(apobec3.transition)'], ax=ax2_1, color='#2E8B57', shade=True, alpha=0.2,zorder=1)
ax1.scatter(np.median(mrcas['age(apobec3.transition)']),167,s=80,color="#2E8B57",zorder=500, edgecolor='black')
ax2_1.set_ylim(0, 7)
ax1.plot([np.median(mrcas['age(apobec3.transition)']), np.median(mrcas['age(apobec3.transition)'])], [0, 167], color='#2E8B57', linestyle='--', linewidth=0.5)

ax1.text(2010, 5, "tMRCA(outgroup)", size=12, color="#0070BB")
ax1.text(2014, 35, "Transition", size=12, color="#2E8B57")
ax1.text(2016.5, 5, "tMRCA(hMPXV-1)", size=12, color="#003262")


ax1.set_xlim(2010, 2025)

for axx in [ax2_1, ax2_2, ax2_5, ax1]:
    axx.spines['left'].set_visible(False)
    axx.spines['right'].set_visible(False)
    axx.spines['top'].set_visible(False)
    axx.yaxis.set_visible(False)


# plt.savefig("Zoom.png", bbox_inches='tight', dpi=700)    


# In[ ]:





