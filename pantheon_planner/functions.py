# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 19:20:08 2018

@author: Alessandro
"""
from __future__ import division
import numpy 

from concorde.tsp import TSPSolver
from numpy import linalg as LA

molt = 1



'''
	my functions!!
'''

def generateCapturePts(trees_pos,n_rows,n_cols,row_step,col_step):
    #assumption: trees are generated row-wise from bottom to up
    
    row_mask = createRowMask2(n_cols,col_step)

    capture_pts = numpy.zeros(shape = ( n_rows * row_mask.shape[0] + n_cols, 2))

    pts = numpy.zeros(shape = ( row_mask.shape[0], 2 ) )

    last_pts = numpy.zeros(shape = ( n_cols, 2 ) )


    for i in range(n_rows):
        curr_row = i*row_step
        pts[:,0] = row_mask
        for j in range(n_cols + 1):
                pts[j,1] = curr_row - row_step/2
  
        for j in range(n_cols,pts.shape[0]):
                pts[j,1] = curr_row
           
        capture_pts[i*row_mask.shape[0]:(i+1)*row_mask.shape[0],:] = pts


    last_row_to_draw = i*row_step
    last_pts[:,0] = row_mask[0:n_cols]

    for j in range(n_cols):
        last_pts[j,1] = last_row_to_draw + row_step/2
    
    capture_pts[-n_cols:,:] = last_pts  

    return capture_pts



def generateinBetweenCapturePts(n_rows,n_cols,row_step,col_step):

    row_mask = createRowMask(n_cols,col_step)

    capture_pts = numpy.zeros(shape = ( (n_rows + 1) * row_mask.shape[0], 2))
    
    pts = numpy.zeros(shape = ( row_mask.shape[0], 2 ) )

    for i in range(n_rows + 1):
        curr_row = i*row_step
        pts[:,0] = row_mask 
        pts[:,1] = curr_row - row_step/2
  
        capture_pts[i*row_mask.shape[0]:(i+1)*row_mask.shape[0],:] = pts

    return capture_pts


def createRowMask(n,step):

    row_mask = numpy.zeros(shape = (n+1))

    row_mask[0] = -step/2
    for i in range(1,n+1,1):
        row_mask[i] = row_mask[i-1] + step
    
    return row_mask

def createRowMask2(n,step):

    row_mask = numpy.zeros(shape = (2*n + 1))

    row_mask[0] = 0

    for i in range(1,n):
        row_mask[i] = row_mask[i-1] + step

    row_mask[n] = -step/2


    for i in range(n+1,row_mask.shape[0]):
        row_mask[i] = row_mask[i-1] + step
    
    return row_mask

def createColMask(n,step):

    col_mask = numpy.zeros(shape = (n+1))

    col_mask[0] = -step/2
    for i in range(1,n+1,1):
        col_mask[i] = col_mask[i-1] + step
    
    return col_mask


def readFile(name,n_nodes):
    #edges = numpy.zeros(shape=(n_nodes, 2),dtype=int)
    nodes = numpy.zeros(shape=(n_nodes, 1),dtype=int)
    file = open(name,"r")
    i=0
    for line in file:
        elem = line.split(" ")
        for e in elem:
            if e!="\n":
                nodes[i,0]=int(e)
                i+=1
         
    '''for i in range(len(nodes)-1):
        edges[i,:] = [nodes[i],nodes[i+1]]
     
    edges[i+1,:] = [nodes[i+1],nodes[0]]'''
    return nodes    
        

def matlab(file):
    ed = numpy.zeros(shape=(400,2),dtype=int)
    
    file = open(file,"r")
    i=0
    for line in file:
        ed[i][0]=line.split(" ")[0]
        ed[i][1]=line.split(" ")[1]
        i+=1
    
    return ed

def getNodes(n):
    nodes = []
    for i in range(n):
        nodes.append(str(i))
    return nodes

def mapPoint(nodes,colors):
    
    dic = {}
    for n in nodes:
        col = colors.get(n)  
        for i in range(4):
            dic[(int(n)*4+i)]=col
            
    return dic       

def generateRandomCenters(num_rows, num_col):
    
    n = num_rows*num_col
    centers = numpy.zeros(shape=(n, 2))
    centers[:,0]=20000 * numpy.random.rand(n)
    centers[:,1]=20000 * numpy.random.rand(n)
    
    
    return centers    

def generateCenters(n_rows,n_col,dist_row,dist_col):
    xx = 0#1*molt
    yy =0#1*molt

    centers = numpy.zeros(shape=(n_rows*n_col, 2))
    h =0;
    for i in range(0,n_rows):  
        for j in range(0,n_col):
            
            centers[h,0] = xx
            centers[h,1] = yy  
            h=h+1
           
            xx =xx +dist_col
            
        yy =yy +dist_row
#        yy = 1*molt
        xx = 0
        
    return centers   

def generateTrees(num_rows, num_cols, dist_row,dist_col):
     
    xx = 0#1*molt
    yy =0#1*molt
    tree_ID = 1

    centers = numpy.zeros(shape=(num_rows*num_cols, 3))
    h =0;
    for i in range(0,num_rows):  
        for j in range(0,num_cols):
            
            centers[h,0] = xx
            centers[h,1] = yy
            centers[h,2] = tree_ID
            h=h+1
            tree_ID = tree_ID + 1
            xx =xx +dist_col
            
        yy =yy +dist_row
#        yy = 1*molt
        xx = 0
        
    return centers        
            
    
def generateRefFrame(centers,r):
    
    h=0
    #1.15176*molt
    print(" Selected tree radious = "+ str(r))
    v = numpy.matrix([[0,r],[r,0],[-r,0],[0,-r]])

    points =numpy.zeros(shape =(4*centers.shape[0],2))
   
    for i in range(centers.shape[0]):
        points[h:(h+4),:]=v+centers[i,:]
        h=h+4

    return points        
            

    
def formatEdges(solution):
    
    edges = numpy.zeros(shape=(solution.size, 2),dtype=int)    
    for i in range(solution.size-1):
        edges[i,:]=[solution[i],solution[i+1]]
        
    edges[i+1,:]=[solution[i+1], solution[0]]    #to close the loop
    return edges
      
def groupByColour(dic):
    dic2 ={}
    
    for key in dic:
        val = dic.get(key)
        if val in dic2:
            elem = dic2[val]
            elem.append(key)
            dic2[val]= elem
        else:
            elem = []
            elem.append(key)
            dic2[val]=elem
    return dic2

def get_Trees_Edges(neighbors,n_trees):
    
    n_edges_i = []
    n_edges_j = []
    h=0
    for key in neighbors:
        elem = neighbors.get(key)
        for i in range(len(elem)):
            n_edges_i.append(int(key))
            n_edges_j.append(int(elem[i]))
            h+=1
    return n_edges_i,n_edges_j


def rotMatrix(x,y,theta,cx,cy):
    p = numpy.zeros(shape=(1, 2))
    new_x = x-cx;
    new_y = y-cy;

    p[0][0]= new_x*numpy.cos(theta)-new_y*numpy.sin(theta)+ cx
    p[0][1]= new_x*numpy.sin(theta)+new_y*numpy.cos(theta)+ cy
    
    return p
    
