# -*- coding: utf-8 -*-
"""
Created on Thu May 23 00:50:09 2013
@author: ryan
"""
import numpy as np

#Function for removing duplicates from an array
def rmv_duplicate(array):
    temp = []
    array.sort()
    for n in range(len(array)-1):
        if(array[n+1] == array[n]):
            temp.append(array[n])
    for m in range(len(temp)):
        array.remove(temp[m])
    return array
    
#Function to remove a subarray from an Array
def rmvNbr1(arr1,arr2):
    for num2 in arr2:
       rmv = []
       for num1 in arr1:
           if(num1 == num2):
               rmv.append(num2)
       for num in rmv:
           arr1.remove(num)
    return arr1

#The Knuth Shuffle Algorithm    
def Knuth(numNodes,array):
    new_arr = []
    limit = numNodes
    for k in range(numNodes):
        index = np.random.randint(limit)
        for n in range(limit):
            if (n == index):
                select_Node = array[index]
                new_arr.append(select_Node)
            elif(n > index):
                m = n-1
                array[m] = array[n]
        limit = limit-1
    return new_arr
    
#Finding the CDS (An implementation of Macker's ECDS Algorithm)
def ECDS(nodes,degree,nbr1,nbr):
       
    #Number of Nodes in the Topology    
    numNodes = len(nodes)    
    '''Start of the ECDS Algorithm'''
    #Initializing the CDS Matrix
    CDS = np.ones(numNodes)
    #Step 1
    #Checking 1 and 2 Hop neighbor of each node
    for m in range(numNodes):
        for n in nbr[m]:
            ptr = nodes.index(n)
            if((degree[m] < degree[ptr]) or (degree[m] == degree[ptr] and m<ptr)):
                CDS[m] = 0
                break
    
    #If degree[m]<2, set CDS[m] = 0
    for m in range(numNodes):
        if(degree[m] < 2):
            CDS[m] = 0
            break
    
    #Step 2 of ECDS Algorithm
    #Part 1: Node with the highest degree in 1-Hop neighborhood
    n1_max = -np.ones(numNodes)
    for m in range(numNodes):
        if(CDS[m] == 0):
            MAX = degree[m]
            for n in nbr1[m]:
                ptr = nodes.index(n)
                if(MAX < degree[ptr]):
                    MAX = degree[ptr]
                    n1_max[m] = n
                elif(MAX == degree[ptr] and m < ptr):
                    n1_max[m] = n
           
    #Part 2: Visiting 1-Hop neighbor of each node and marking them
    #Finding if route exist from 1-Hop neighbor with Max Degree to all
    #other 1-Hop neighbor of node "n" wihtout traversing an intermediate
    #node of degree less than "degree[n]"
    
    #Visited Matrix of 1 and 2-Neighbors
    visitor = np.zeros([numNodes, numNodes])
    for m in range(numNodes):
        for n in nbr[m]:
            ptr = nodes.index(n)
            visitor[m][ptr] = 1 
    
    for m in range(numNodes):
        #Initializing the Stack
        stack = []
        if(CDS[m] == 0) and (n1_max[m] != -1):
            stack.append(n1_max[m])
            while(len(stack) > 0):
                K = int(stack.pop())
                k = nodes.index(K)
                visitor[m][k] = 0            
                for n in nbr1[k]:
                    ptr = nodes.index(n)
                    if(visitor[m][ptr] == 1 and ptr != m):
                        visitor[m][ptr] = 0
                        if(degree[m] < degree[ptr] or (degree[m] == degree[ptr] and m<ptr)):
                            stack.append(n)
                                                    
    for m in range(numNodes):
        if(CDS[m] == 0):
            for n in nbr1[m]:
                ptr = nodes.index(n)
                if(visitor[m][ptr] == 1):
                    CDS[m] = 1
        
    return CDS
'''End of the ECDS Algorithm'''

#Dijkstra's Shortest Path Algorithm
def Dijkstra(source,num_Nodes,bigLamda,Nodes,wt_color,Nbr1):
    wt_Lamda = np.ones(num_Nodes)
    temp = bigLamda.min()
    if(temp != 0):
        weight2 = bigLamda/temp
        for n in Nodes:
            ptr = Nodes.index(n)
            wt_Lamda[n] = weight2[ptr]    
        
    #Holds the nodes which forms the shortest path
    prev = -np.ones(num_Nodes)
    #The array of visitor nodes
    V = np.zeros(num_Nodes)
            
    #Array to hold the distance of a node from the source
    dis = np.ones(num_Nodes)*num_Nodes*2
    ptr = Nodes.index(source)
    for node in Nbr1[ptr]:
        dis[node] = 1*wt_color[node]
            
    #The source is marked as visited
    V[source] = 1
    Num_Nodes = len(Nodes)          
    for count in range(Num_Nodes):
        index = -2
        min_distance = num_Nodes
        for node in Nodes:
            if(V[node] == 0 and dis[node] < min_distance):
                min_distance = dis[node]
                index = node                            
        if(index != -2):        
            V[index] = 1
            ptr = Nodes.index(index)
            for node in Nbr1[ptr]:
                dist = 1*wt_color[node]
                if(dis[index] + dist < dis[node]):
                    dis[node] = dis[index] + dist
                    prev[node] = index                        
    return prev
    
#The min hop path for a node
def minHopPath(dest,prev,array):
    if(prev[dest] != -1):
        minHopPath(prev[dest],prev,array)   
    array.append(dest)
    return array
    
#The Associated Source Nodes
def SourceNodeArray(num_Nodes,Source,Load,CDS,nbr1,Nodes,max_client):
    for n in range(num_Nodes):
        if(CDS[n] == 0 and len(nbr1[n]) >=1):
            temp_arr = []
            load = max_client
            for n1 in nbr1[n]:
                if(CDS[n1] == 1):
                    temp_arr.append(n1)
            if(len(temp_arr) > 0):        
                for t1 in temp_arr:
                    ptr = Nodes.index(t1)
                    if(load >= Load[ptr]):    
                        load = Load[ptr]
                        source = t1
                ptr = Nodes.index(source)
                Load[ptr] += 1
                Source[n] = source                    
        elif (CDS[n] == 1):
            Source[n] = n
    for n in range(num_Nodes):
        temp_arr = []
        if(Source[n] == -1):
            for n1 in nbr1[n]:
                if(CDS[n1] == 1):
                    Source[n] = n1
                        
    return Source