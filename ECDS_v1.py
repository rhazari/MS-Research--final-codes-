# -*- coding: utf-8 -*-
"""
Created on Sun Sep 01 15:56:45 2013
@author: ryan
"""

import numpy as np
import pylab as pl
import math
import customFuncs as cf
import time
#import sys
#import pdb

start = time.time()
#np.random.seed(62)
#np.random.seed(55)

#Number of nodes
num_Nodes = 200
#num_Nodes = 25

#Communication Range:The 1-Hop Neighbor Range
#nbr_range = 380
nbr_range = 300


#Max number of clients per dominant nodes
max_client = 10

#Length and Breadth of the Square Area
length = 10**3; breadth = 10**3

#X and Y cordinates of the nodes in the square area
x_cor = np.random.rand(num_Nodes)*length; y_cor = np.random.rand(num_Nodes)*breadth

#The Nodes
nodes=[]
for m in range(num_Nodes):
    nodes.append(m)
    pl.plot(x_cor[m],y_cor[m],'b+')

#Array holding X & Y cordinates of each node
nodeArray = []
for i in range(num_Nodes):
    nodeArray.append([x_cor[i],y_cor[i]])
    
#Function for finding distance between 2 cordinates    
def distance(node1,node2):
    x1 = nodeArray[node1][0]; x2 = nodeArray[node2][0]
    y1 = nodeArray[node1][1]; y2 = nodeArray[node2][1]
    path_length = math.sqrt((x1-x2)**2 + (y1-y2)**2)
    if (path_length < nbr_range):
        return path_length
    else:
        return length

#Function to Find 1-Hop neighbors
def OneHopNbr(arr1):
    arr2 = []
    for m in arr1:
        subarr = []
        for n in arr1:
            if(m != n and distance(m,n) < nbr_range):
                subarr.append(n)
        arr2.append(subarr)
    return arr2

#Function to find 2-Hop Neighbors 
def TwoHopNbr(arr1,arr2):
    arr3 = []
    for m in arr1:
        subarr = []
        ptr1 = arr1.index(m)
        for n1 in arr2[ptr1]:
            ptr2 = arr1.index(n1)
            for n2 in arr2[ptr2]:
                if(m != n2):
                    subarr.append(n2)
            cf.rmv_duplicate(subarr)
            cf.rmvNbr1(subarr,arr2[ptr1])
        arr3.append(subarr)
    return arr3

#Array of 1-Neighbors of all the nodes
nbr1 = []
nbr1 = OneHopNbr(nodes)

#Degree of each Node(i.e the no. of 1 Hop neighbors)
degree = np.zeros(num_Nodes)
for node in nodes:
    degree[node] = len(nbr1[node])

#Array 2-Neighbors of all the nodes 
nbr2 = []
nbr2 = TwoHopNbr(nodes,nbr1)

#Array of 1 and 2-Hop Neighbors of all the nodes
nbr = []
for m in range(num_Nodes):
    nbr.append(nbr1[m]+nbr2[m])

#Array holder of the dominant nodes    
CDS = cf.ECDS(nodes,degree,nbr1,nbr)
    
#Array of nodes in the CDS
Nodes = []
Num_Nodes = 0
for m in range(num_Nodes):
    if(CDS[m]==1):
        Num_Nodes += 1
        Nodes.append(m)
           
#Array of 1-Neighbors of nodes in the CDS
Nbr1 = []
Nbr1 = OneHopNbr(Nodes)
    
#Array 2-Neighbors of nodes in the CDS 
Nbr2 = []
Nbr2 = TwoHopNbr(Nodes,Nbr1)
    
#Array of 1 and 2-Neighbors of nodes in the CDS
Nbr = []
for m in range(Num_Nodes):
    Nbr.append(Nbr1[m]+Nbr2[m])
        
#Calling on the Knuth Shuffle
FY_Nodes = cf.Knuth(Num_Nodes,Nodes)
    
#Array of nodes in the CDS(after the Knuth Shuffle)
Nodes = []
for m in range(num_Nodes):
    if(CDS[m] == 1):
        Nodes.append(m)
    
colorArray = np.zeros(Num_Nodes)  
#Function to color a node
def colorNode(node):
    subarr = []; colors = []
    ptr1 = Nodes.index(node)
    for c in range(Num_Nodes):
        colors.append(c+1)
    for n1 in Nbr[ptr1]:
        ptr2 = Nodes.index(n1)
        if(colorArray[ptr2] != 0):
            subarr.append(colorArray[ptr2])
    cf.rmvNbr1(colors,subarr)
    colors.reverse()
    colorArray[ptr1] = colors.pop()
        
#Coloring Nodes in CDS after the Knuth Shuffle    
for m in range(Num_Nodes):
    fy_node = FY_Nodes[m]
    colorNode(fy_node)

#Finding Maximum Color number in the 1 & 2 Neighborhood of each node
maxColor = []
for n1 in range(Num_Nodes):
    color = colorArray[n1]
    for n2 in Nbr1[n1]:
        ptr = Nodes.index(n2)
        if(color < colorArray[ptr]):
            color = colorArray[ptr]
    maxColor.append(color)
        
#Total Transmission Slot for each node
FrameSize = []
for n1 in range(Num_Nodes):
    color = maxColor[n1]
    power = math.ceil(math.log(color)/math.log(2))
    FrameSize.append(2**power)
    
#Calculating the maximum frame size
maxColorNum = colorArray.max()
power = math.ceil(math.log(maxColorNum)/math.log(2))
maxFrame = int(2**power)
    
#Matrix for Tx slot of each node (based on color numbers) 
TxSlot = np.zeros([Num_Nodes,maxFrame])
                    
#Matrix for Tx slot of nodes
Tx = np.zeros([Num_Nodes,maxFrame])
                
#Initializing the TxSlot Matrix    
for n1 in range(Num_Nodes):
    color = colorArray[n1]
    power = math.ceil(math.log(color)/math.log(2)) 
    n = 0; slot = color
    while(slot <= maxFrame):
        n += 1
        TxSlot[n1][slot-1] = 1
        slot = color + n*(2**power)    
                
#Tx Slot for the Dominant Nodes     
for n1 in range(Num_Nodes):
    Tx[n1] = TxSlot[n1]
    for n2 in Nbr[n1]:
        ptr = Nodes.index(n2)
        if(colorArray[n1] < colorArray[ptr]):
            Tx[n1] = Tx[n1]-TxSlot[ptr]
                
#Calculating the number of Tx slots in a Frame for each Dominant Node 
TxCount = np.zeros(Num_Nodes)
for n1 in range(Num_Nodes):
    for n2 in range(maxFrame):
        if(Tx[n1][n2] == 1):
            TxCount[n1] += 1
                
#Transmission % for each of the dominant nodes
TxPct = TxCount/maxFrame 

#Weights for dominant nodes (inverse of the Transmission %) 
invTxPct = maxFrame/TxCount
invTxPct = invTxPct/invTxPct.min()            
wt_color = np.zeros(num_Nodes)
for n in Nodes:
    ptr = Nodes.index(n)
    wt_color[n] = invTxPct[ptr]
 
#Load Count...
Load = np.ones(Num_Nodes)                           
#The Associated Source Nodes
Source = -np.ones(num_Nodes)
cf.SourceNodeArray(num_Nodes,Source,Load,CDS,nbr1,Nodes,max_client)
                           
#Initializing the BigLamda
#The path returned by Dijkstra does not include the Source Node
#Each node generates traffic for every other node in the network
bigLamda = np.zeros(Num_Nodes)
                   
#The Main
#Calculating the bigLamda for the nodes in CDS                   
for n in range(num_Nodes):
    prevPath = []
    source = int(Source[n])
    if(source != -1):
        #Dijkstra shortest path algorithm...
        prevPath = cf.Dijkstra(source,num_Nodes,bigLamda,Nodes,wt_color,Nbr1)
        #The destination nodes (computing the BigLamda)                     
        for m in range(num_Nodes):
            path = []
            dest = -1
            if(Source[m] != source):
                dest = Source[m]                                            
            if(dest != -1):   
                cf.minHopPath(dest,prevPath,path)               
                path.pop()
                path.append(source)
                num_hops = len(path)
                if(num_hops >= 1):
                    for n3 in range(num_hops):
                        node = path.pop()
                        ptr = Nodes.index(node)
                        bigLamda[ptr] = bigLamda[ptr]+1
                
#The LoadFactor...            
LoadFactor = np.zeros(Num_Nodes);
for n1 in range(Num_Nodes):
    LoadFactor[n1] = bigLamda[n1]*maxFrame/TxCount[n1]
    
#Writing Data into File...    
'''my_file = open(sys.argv[1],'w')
my_file.write(str(float(Num_Nodes)/num_Nodes)+"\n")
my_file.write(str(LoadFactor.max())+"\n")
my_file.write(str(LoadFactor.min())+"\n")
my_file.write(str(LoadFactor.sum()/Num_Nodes)+"\n")
my_file.close()'''
    
#Print Details
print "Fraction of nodes in the CDS backbone = "+str(float(Num_Nodes)/num_Nodes)       
print "Max LoadFactor :"+str(LoadFactor.max())
print "Min LoadFactor :"+str(LoadFactor.min())
print "Avg LoadFactor :"+str(LoadFactor.sum()/Num_Nodes)
print Load

sel_LF = []; sel_bigLamda = []; sel_TxPct = []; sel_color = []; sel_Cor = [];

maxLF = LoadFactor.max()
for k in range(Num_Nodes):
    if(LoadFactor[k] == maxLF):
        ptr1 = k
        
for node in Nbr[ptr1]:
    ptr2 = Nodes.index(node)
    sel_LF.append(LoadFactor[ptr2])
    sel_bigLamda.append(bigLamda[ptr2])
    sel_TxPct.append(TxPct[ptr2])
    sel_color.append(colorArray[ptr2])
    sel_Cor.append(nodeArray[node])
    #pl.plot(x_cor[node],y_cor[node],'ro')

end = time.time()
print end - start

#Plot and Display
for node in Nodes:
    pl.plot(x_cor[node],y_cor[node],'go')

pl.plot(x_cor[Nodes[ptr1]],y_cor[Nodes[ptr1]],'rs')
            
for m in Nodes:
    ptr = Nodes.index(m)
    for n in Nbr1[ptr]:
        pl.plot([x_cor[m], x_cor[n]],[y_cor[m], y_cor[n]],'k:')    
                                    
pl.xlim(-0.2,length+0.2)
pl.ylim(-0.2,length+0.2)
pl.show()
#pl.savefig('ECDS.pdf', dpi=1000)'''