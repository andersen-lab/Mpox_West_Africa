#!/usr/bin/env python
# coding: utf-8

# In[121]:


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


# In[146]:


treeFile="hmpxv.nex" 
ll=bt.loadNexus(treeFile)
ll.traverse_tree()
ll.treeStats()



ancestor=ll.commonAncestor(ll.getExternal(lambda k: k.name in 
                                          ['unpub|TRM320|Nigeria|Rivers|2022-10-03', 
                                           'unpub|TRM289|Nigeria|Abia|2022-09-25'])) 
ll.collapseSubtree(ancestor,'A.2.3',widthFunction=lambda x:3)

ancestor2=ll.commonAncestor(ll.getExternal(lambda k: k.name in 
                                          ['OP535316|MPXV_5218|Nigeria|2018-02', 
                                           'MK783029|3029|Nigeria|2017-12-06'])) 
ll.collapseSubtree(ancestor2,'A.2.1',widthFunction=lambda x:3)

ancestor3=ll.commonAncestor(ll.getExternal(lambda k: k.name in 
                                          ['OP612676|MPXV/013/19|Nigeria|2019-02-26', 
                                           'OP535337|MPXV_2945|Nigeria|2017-10'])) 
ll.collapseSubtree(ancestor3,'A.2.2',widthFunction=lambda x:3)

ancestor4=ll.commonAncestor(ll.getExternal(lambda k: k.name in 
                                          ['MK783033|2920|Nigeria|2017-10-09', 
                                           'unpub|TRM283|Nigeria|Ebonyi|2022'])) 
ll.collapseSubtree(ancestor4,'A.2.4',widthFunction=lambda x:3)

ancestor5=ll.commonAncestor(ll.getExternal(lambda k: k.name in 
                                          ['MT903343|MPXV-UK_P1|UK_ex_Nigeria|2018-09-07', 
                                           'ON563414|MA001|USA|2022-05-19'])) 
ll.collapseSubtree(ancestor5,'A.2.5',widthFunction=lambda x:3)



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
        
for xx in range(1955, 2025, 10):
    plt.axvspan(xx, xx+5, facecolor='#F5F5F5', alpha=1, zorder=1, edgecolor=None)  



ll.plotTree(ax1,x_attr=x_attr,colour=c_func,zorder=100, width=1.2) 
ll.plotPoints(ax1,x_attr=x_attr,size=35,colour=c_func2,zorder=100) 


for k in ll.Objects: 
    if isinstance(k,bt.clade): 
        x=2023.4137-k.traits['height']
        y=k.y
        clade=plt.Polygon(([x,y-0.001*len(ll.Objects)],
                           [x,y+0.001*len(ll.Objects)],
                           [2024,y+k.width/2.0],
                           [2024,y-k.width/2.0]),
                          facecolor='#89CFF0',edgecolor="black",zorder=1) 

        ax1.add_patch(clade)


ax2_1 = ax1.twinx()
sns.kdeplot(mrcas['tMRCA(1971)'], ax=ax2_1, color='#8B0000', shade=True, alpha=0.2,zorder=1)
ax1.scatter(np.median(mrcas['tMRCA(1971)'])-0.1,30,s=30,color="#8B0000",zorder=500, edgecolor='black')
ax2_1.set_ylim(0, 1.5)
ax1.plot([np.median(mrcas['tMRCA(1971)'])-0.1, np.median(mrcas['tMRCA(1971)'])-0.1], [0, 30], color='#8B0000', linestyle='--', linewidth=0.5)

ax2_3 = ax1.twinx()
sns.kdeplot(mrcas['tMRCA(akwa)'], ax=ax2_3, color='#ED9121', shade=True, alpha=0.2,zorder=1)
ax1.scatter(np.median(mrcas['tMRCA(akwa)'])-0.1,38.8,s=30,color="#ED9121",zorder=500, edgecolor='black')
ax2_3.set_ylim(0, 1)
ax1.plot([np.median(mrcas['tMRCA(akwa)'])-0.1, np.median(mrcas['tMRCA(akwa)'])-0.1], [0, 38.8], color='#ED9121', linestyle='--', linewidth=0.5)


ax2_4 = ax1.twinx()
sns.kdeplot(mrcas['tMRCA(Cameroon)'], ax=ax2_4, color='#C04000', shade=True, alpha=0.2,zorder=1)
ax1.scatter(np.median(mrcas['tMRCA(Cameroon)'])-0.1,40.8,s=30,color="#C04000",zorder=500, edgecolor='black')
ax2_4.set_ylim(0, 1)
ax1.plot([np.median(mrcas['tMRCA(Cameroon)'])-0.1, np.median(mrcas['tMRCA(Cameroon)'])-0.1], [0, 40.8], color='#C04000', linestyle='--', linewidth=0.5)


ax1.set_xlim(1960, 2025)


ax1.text(1972.5, 32.05, "Abia 1971", size=10)
ax1.text(2024.5, 30.5, "Zx", size=12)

ax1.scatter(1972,20,s=200,color='#FFA000',zorder=1, edgecolor='black')
ax1.text(1974, 19.5, "Cameroon", size=12)

ax1.scatter(1972,18,s=200,color='#00356B',zorder=1, edgecolor='black')
ax1.text(1974, 17.5, "Abia", size=12)

ax1.scatter(1972,16,s=200,color='#4997D0',zorder=1, edgecolor='black')
ax1.text(1974, 15.5, "Akwa Ibom", size=12)


ax1.scatter(1972,14,s=200,color='#89CFF0',zorder=1, edgecolor='black')
ax1.text(1974, 13.5, "hMPXV-1", size=12)


for axx in [ax2_1, ax2_2, ax2_3, ax2_4, ax2_5, ax1]:
    axx.spines['left'].set_visible(False)
    axx.spines['right'].set_visible(False)
    axx.spines['top'].set_visible(False)
    axx.yaxis.set_visible(False)



ax1.text(1970, 5, "tMRCA(KJ642617, hMPXV-1)", size=7, color="#8B0000")
ax1.text(1989, 3, "tMRCA(Cameroon)", size=7, color="#C04000")
ax1.text(1994, 8, "tMRCA(Awka Ibom, Mbongue)", size=7, color="#ED9121")


clade=plt.Polygon(([2013,-1],
                   [2013,32],
                   [2025,32],
                   [2025, -1]),

                  facecolor='#3AA8C1',edgecolor=None,zorder=10000, alpha=0.1) 

ax1.add_patch(clade)



# plt.savefig("Full.png", bbox_inches='tight', dpi=700)    

