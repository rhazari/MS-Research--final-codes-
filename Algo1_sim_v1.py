# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:28:51 2013
@author: ryan
"""

import sys
import numpy as np
import math
import customFuncs as cf
#import pdb

#Number of nodes
num_Nodes = 100

#Max Contention Load
max_client = 10

#Length and Breadth of the Square Area
#length = math.sqrt(num_Nodes); breadth = math.sqrt(num_Nodes)
length = 10; breadth = 10

#Communication Range:The 1-Hop neighbor Range
OneHopRange = np.linspace(2.0,4.5,15)
count = len(OneHopRange)

#Parameter to find
MaxLoadFactor = np.zeros(count)
MinLoadFactor = np.zeros(count)
AvgLoadFactor = np.zeros(count)
#percent_node_contention = np.zeros(count)
#Sel_bigLamda = np.zeros(count)
#Avg_bigLamda = np.zeros(count)
#Sel_TxPct = np.zeros(count)
#Avg_TxPct = np.zeros(count)
Num_CDS_Nodes = np.zeros(count)

#Number of Test Runs for each 1Hop Neighbor Range
test_runs = 5

for a in range(count):
    nbr_range = OneHopRange[a]
    
    for b in range(test_runs):
        #X and Y cordinates of the nodes in the square area
        x_cor = np.random.rand(num_Nodes)*length; y_cor = np.random.rand(num_Nodes)*breadth
        
        #The Nodes
        nodes=[]
        for m in range(num_Nodes):
            nodes.append(m)
        
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
                return num_Nodes
        
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
        read = 1; sel_node = -1; 
        maxLF = 0; color_flag = 0 
        sel_bck_node = -1; 
        counter_1 = 0; counter_2 = 0
        bck_node_limit = num_Nodes; 
        LF_val_limit = 3
        while(read):
            if(sel_node == -1):        
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
                if(color_flag == 0):
                    color_flag = 1
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
                    colorNode(FY_Nodes[m])
                    
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
                            
                #Computing the Network Diameter    
                maxHop = 0        
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
                                if(maxHop < len(path)):
                                    maxHop = len(path)                   
                            if(dest != -1):
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
                
                #Dominant Nodes with their LoadFactors...
                LF_Nodes = []    
                for n1 in range(Num_Nodes):
                    LF_Nodes.append([LoadFactor[n1],Nodes[n1]])    
                #Sorting in the reverse based on LoadFactors of the Dominant Nodes
                LF_Nodes.sort() 
                LF_Nodes.reverse()
            
            if(maxLF == 0):
                maxLF = LoadFactor.max()
            if(sel_node == -1 and counter_2 == 0):
                sel_node = LF_Nodes[counter_2][1]
            
            #Storing the max LoadFactors over the different iterations of 1-Hop Neighbor(non-Dominant Nodes)
            #Iteratively finding the non-Dominant Node which minimizes the MaxLoadFactor
            maxLF_arr = []
            for node in nbr1[sel_node]:
                if(CDS[node] != 1 and counter_1 < bck_node_limit):
                    
                    Nodes_t = []; Num_Nodes_t = Num_Nodes+1
                    Color_Array = np.zeros(Num_Nodes_t)
                    for m in range(Num_Nodes):
                        Color_Array[m] = colorArray[m]
                        Nodes_t.append(Nodes[m])
                    
                    Nodes_t.append(node)
                    Nbr1_t = []; Nbr1_t = OneHopNbr(Nodes_t)
                    Nbr2_t = []; Nbr2_t = TwoHopNbr(Nodes_t,Nbr1_t)
                    
                    Nbr_t = []
                    for m in range(Num_Nodes_t):
                        Nbr_t.append(Nbr1_t[m]+Nbr2_t[m])
                    
                    CDS_t = np.zeros(num_Nodes)
                    for n1 in Nodes_t:
                        CDS_t[n1] = 1        
                    
                    #Coloring of the Node...
                    color = []; subarr = []
                    for c in range(Num_Nodes_t):
                        color.append(c+1)
                    ptr1 = Nodes_t.index(node)
                    for n in Nbr_t[ptr1]:
                        ptr2 = Nodes_t.index(n)
                        if(colorArray[ptr2] != 0):
                            subarr.append(colorArray[ptr2])  
                    cf.rmvNbr1(color,subarr)
                    color.reverse()
                    Color_Array[Num_Nodes] = color.pop()
                    
                    #Calculating the maximum frame size
                    maxColor_t = Color_Array.max()
                    power = math.ceil(math.log(maxColor_t)/math.log(2))
                    maxFrame_t = int(2**power)
                    
                    #Matrix for slot Tx of each node
                    TxSlot_t = np.zeros([Num_Nodes_t,maxFrame_t])
                
                    #Matrix for slot Tx of nodes
                    Tx_t = np.zeros([Num_Nodes_t,maxFrame_t])
                    
                    #Initializing the TxSlot Matrix    
                    for n1 in range(Num_Nodes_t):
                        color = Color_Array[n1]
                        power = math.ceil(math.log(color)/math.log(2)) 
                        n = 0; slot = color
                        while(slot <= maxFrame_t):
                            n += 1
                            TxSlot_t[n1][slot-1] = 1
                            slot = color + n*(2**power)
                    
                    #Tx Slot for the Dominant Nodes    
                    for n1 in range(Num_Nodes_t):
                        Tx_t[n1] = TxSlot_t[n1]
                        for n2 in Nbr_t[n1]:
                            ptr = Nodes_t.index(n2)
                            if(Color_Array[n1] < Color_Array[ptr]):
                                Tx_t[n1] = Tx_t[n1]-TxSlot_t[ptr]
                    
                    #Calculating the number of Tx slots in a Frame for each Dominant Node 
                    TxCount_t = np.zeros(Num_Nodes_t)
                    for n1 in range(Num_Nodes_t):
                        for n2 in range(maxFrame_t):
                            if(Tx_t[n1][n2] == 1):
                                TxCount_t[n1] += 1
                    
                    #Weights for dominant nodes (inverse of the Transmission %)            
                    invTxPct_t = maxFrame_t/TxCount_t
                    invTxPct_t = invTxPct_t/invTxPct_t.min()                    
                    wt_clr = np.zeros(num_Nodes)
                    for n in Nodes_t:
                        ptr = Nodes_t.index(n)
                        wt_clr[n] = invTxPct_t[ptr]
                        
                    #Load Count...
                    Load_t = np.ones(Num_Nodes_t)                                    
                    #The Associated Source Nodes
                    Source_t = -np.ones(num_Nodes)
                    cf.SourceNodeArray(num_Nodes,Source_t,Load_t,CDS_t,nbr1,Nodes_t,max_client)
                 
                    #Initializing the BigLamda
                    #The path returned by Dijkstra does not include the Source Node
                    #Each node generates traffic for every other node in the network
                    bigLamda_t = np.zeros(Num_Nodes_t)
                                      
                    #The Main
                    #Calculating the bigLamda for the nodes in CDS                   
                    for n in range(num_Nodes):
                        prevPath = []
                        source = int(Source_t[n])
                        if(source != -1):
                            prevPath = cf.Dijkstra(source,num_Nodes,bigLamda_t,Nodes_t,wt_clr,Nbr1_t)                               
                            for m in range(num_Nodes):
                                path = []
                                dest = -1
                                if(Source_t[m] != source):
                                    dest = Source_t[m]        
                                        
                                if(dest != -1):   
                                    cf.minHopPath(dest,prevPath,path)                                       
                                if(dest != -1):
                                    path.pop()
                                    path.append(source)
                                    num_hops = len(path)
                                    if(num_hops >= 1):
                                        for n1 in range(num_hops):
                                            n2 = path.pop()
                                            ptr = Nodes_t.index(n2)
                                            bigLamda_t[ptr] = bigLamda_t[ptr]+1
                                    
                    #The LoadFactor...            
                    LoadFactor_t = np.zeros(Num_Nodes_t);
                    for n1 in range(Num_Nodes_t):
                        LoadFactor_t[n1] = bigLamda_t[n1]*maxFrame_t/TxCount_t[n1]
                    
                    maxLF_t = LoadFactor_t.max()
                    maxLF_arr.append([maxLF_t,node])
            
            if(counter_1 < bck_node_limit):        
                new_maxLF = maxLF
                for item in maxLF_arr:
                    if(item[0] < new_maxLF):
                        new_maxLF = item[0]
                        sel_bck_node = item[1]
                           
                if(new_maxLF == maxLF and counter_2 < LF_val_limit):
                    counter_2 += 1
                    sel_node = LF_Nodes[counter_2][1]
                    sel_bck_node = -1
                elif(new_maxLF < maxLF and counter_2 <= LF_val_limit):
                    maxLF = new_maxLF
                    CDS[sel_bck_node] = 1
                    FY_Nodes.append(sel_bck_node)
                    counter_1 += 1
                    sel_node = -1; counter_2 = 0
                elif(new_maxLF == maxLF and counter_2 == LF_val_limit):
                    read = 0            
            else:
                read = 0
             
        #Collecting Statistics
        #The node with the Max LoadFactor
        maxLF = LoadFactor.max()
        for k in range(Num_Nodes):
            if(LoadFactor[k] == maxLF):
                ptr1 = k
                
        MaxLoadFactor[a] += LoadFactor.max()
        MinLoadFactor[a] += LoadFactor.min()
        AvgLoadFactor[a] += LoadFactor.sum()/Num_Nodes
        #Sel_bigLamda[a] += bigLamda[ptr1]
        #Avg_bigLamda[a] += bigLamda.sum()/Num_Nodes
        #Sel_TxPct[a] += TxPct[ptr1]
        #Avg_TxPct[a] += TxPct.sum()/Num_Nodes
        Num_CDS_Nodes[a] += Num_Nodes
        
        '''contention_nodes = 0
        for m in range(Num_Nodes):
            if(Load[m] >= max_client):
                contention_nodes += 1        
        percent_node_contention[a] += float(contention_nodes)/Num_Nodes'''        
        
MaxLoadFactor = MaxLoadFactor/test_runs
MinLoadFactor = MinLoadFactor/test_runs
AvgLoadFactor = AvgLoadFactor/test_runs
#percent_node_contention = percent_node_contention/test_runs
#Sel_bigLamda = Sel_bigLamda/test_runs
#Sel_TxPct = Sel_TxPct/test_runs
Num_CDS_Nodes = Num_CDS_Nodes/test_runs

my_file = open(sys.argv[1],'w')
for m in range(count):
    my_file.write(str(MaxLoadFactor[m])+"\t")
    
for m in range(count):
    my_file.write(str(MinLoadFactor[m])+"\t")
    
for m in range(count):
    my_file.write(str(AvgLoadFactor[m])+"\t")

for m in range(count):
    my_file.write(str(Num_CDS_Nodes[m])+"\t") 
my_file.close()

'''for m in range(count):
    my_file.write(str(avg_MaxLF[m])+"\t")
    
for m in range(count):
    my_file.write(str(Sel_bigLamda[m])+"\t")
    
for m in range(count):
    my_file.write(str(Sel_TxPct[m])+"\t")'''    


