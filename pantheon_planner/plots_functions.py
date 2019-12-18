#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:23:33 2019

@author: Jacopo
"""

import matplotlib.pyplot as plt

def plotMapTrees(trees_pos,tree_radious):
    
    plt.figure(figsize=(8,6))
    plt.scatter(trees_pos[:,0],trees_pos[:,1],color='g',s=100 * tree_radious)
    
    for i, txt in enumerate(trees_pos[:,2]):
        plt.annotate(txt, (trees_pos[i,0],trees_pos[i,1]))
    
    plt.axis('equal')
    plt.title("Trees positions")
    plt.show()

def plotCapturePts(capture_pos):
    
    plt.figure(figsize=(8,6))
    plt.scatter(capture_pos[:,0],capture_pos[:,1],color='b')
    
    plt.axis('equal')
    plt.title("Captures points positions")
    plt.show()


def plotTreesAndCaptures(trees_pos,capture_pos,tree_radious):
    
    plt.figure(figsize=(8,6))
    ax = plt.gca()
    ax.cla() 
#    plt.scatter(trees_pos[:,0],trees_pos[:,1],color='g', s = 100 * tree_radious)
    plt.scatter(capture_pos[:,0],capture_pos[:,1],color='b')

    for i, txt in enumerate(trees_pos[:,2]):
        circle3 = plt.Circle((trees_pos[i,0], trees_pos[i,1]), tree_radious, color='g', clip_on=False)
        plt.annotate(int(txt), (trees_pos[i,0],trees_pos[i,1]))
        ax.add_artist(circle3)

    for i, txt in enumerate(capture_pos[:,0]):
        plt.annotate(str(i+1), (capture_pos[i,0],capture_pos[i,1]))
    
    plt.axis('equal')
    plt.title("All the points")

    name = "capture_pts_and_trees"
    plt.savefig("/home/majo/libs/pyconcorde/pantheon_planner/results/"+ name, dpi = 800 ,format = 'pdf')
    plt.show()


def plotColouredGraph(centers,neighbors,colors,n_trees):
    x = centers[:,0].flatten()
    y = centers[:,1].flatten()
    
    dic = groupByColour(colors)
    reds = dic.get('R')
    greens = dic.get('G') 
    blues =dic.get('B')
    #yellows = dic.get('Y')
    
    plt.figure(figsize=(8,6))
    for i in range(len(reds)):
        plt.plot(centers[int(reds[i]),0],centers[int(reds[i]),1],'o',color ='r',markersize=15)
        
    for i in range(len(greens)):
        plt.plot(centers[int(greens[i]),0],centers[int(greens[i]),1],'o',color ='g',markersize=15)
        
    for i in range(len(blues)):
        plt.plot(centers[int(blues[i]),0],centers[int(blues[i]),1],'o',color ='b',markersize=15)
    
    '''for i in range(len(yellows)):
        plt.plot(centers[int(yellows[i]),0],centers[int(yellows[i]),1],'o',color ='y',markersize=2)'''
    
    e1,e2 = get_Trees_Edges(neighbors,n_trees)
    e_x = [x[e1],x[e2]]
    e_y = [y[e1],y[e2]]
    plt.plot(e_x, e_y, linestyle='-', color='y')  
 




def plotPoints(points,edges):
    x = points[:,0].flatten()
    y = points[:,1].flatten()
    fig, ax= plt.subplots(figsize=(10,8))
    #ax.scatter(points[:,0],points[:,1])
    
    num =range(points.shape[0])  
    
    for i,txt in enumerate(num):
        ax.annotate(txt, (points[i,0], points[i,1]),)
    
    plt.plot(x[edges.T], y[edges.T], linestyle='-', color='y')     

def plotWithoutLines(points,centers,n_col,n_row):
    x = points[:,0].flatten()
    y = points[:,1].flatten()
    plt.figure(figsize=(16,12))
    plt.plot(centers[:,0],centers[:,1],'o',color ='g')
    plt.plot(x,y,'o',color ='r')
    name = "points "+str(n_col)+"x"+str(n_row)
    plt.savefig("/home/majo/libs/pyconcorde/pantheon_planner/results/"+name, dpi = 800 ,format = 'pdf')
    

def plotSolutionGraph(capture_pos,edges,trees_pos,filename):
    
    x = capture_pos[:,0]
    y = capture_pos[:,1]

    plt.figure(figsize=(8,6))
    ax = plt.gca()
    ax.cla() 

    print(edges)

#    plt.figure(figsize=(12,8))
#    plt.plot(trees_pos[:,0],trees_pos[:,1],'o',color ='g')

    plt.scatter(capture_pos[:,0],capture_pos[:,1],color='b')

    for i, txt in enumerate(trees_pos[:,2]):
        circles = plt.Circle((trees_pos[i,0], trees_pos[i,1]), 0.3, color='g', clip_on=False)
        plt.annotate(int(txt), (trees_pos[i,0],trees_pos[i,1]))
        ax.add_artist(circles)

    for i, txt in enumerate(capture_pos[:,0]):
        plt.annotate(str(i+1), (capture_pos[i,0],capture_pos[i,1]))

    plt.plot(x[edges.T], y[edges.T], linestyle='-', color='r') 

    plt.annotate("start/stop", (x[edges[0,0]]-0.5,y[edges[0,0]]-0.5))

    plt.savefig("/home/majo/libs/pyconcorde/pantheon_planner/results/"+filename, format = 'pdf')

    plt.axis("equal")
    plt.show()


def plotGraph(points,edges,trees_pos,n_col,n_row):
    
    x = points[:,0].flatten()
    y = points[:,1].flatten()

    plt.figure(figsize=(8,6))
    ax = plt.gca()
    ax.cla() 

#    plt.figure(figsize=(12,8))
#    plt.plot(trees_pos[:,0],trees_pos[:,1],'o',color ='g')

    for i, txt in enumerate(trees_pos[:,2]):
        circles = plt.Circle((trees_pos[i,0], trees_pos[i,1]), 0.3, color='g', clip_on=False)
        plt.annotate(txt, (trees_pos[i,0],trees_pos[i,1]))
        ax.add_artist(circles)


    plt.plot(x[edges.T], y[edges.T], linestyle='-', color='y',markerfacecolor='red', marker='o') 
    plt.plot(x,y,'o',color ='r')
    plt.plot(x[edges.T], y[edges.T], color='y',markerfacecolor='red', marker='o') 
    
    name = str(n_col)+"x"+str(n_row)
    plt.savefig("/home/majo/libs/pyconcorde/pantheon_planner/results/"+name, format = 'pdf')
    

    plt.axis("equal")

    plt.show()
    
def plotGraphLight(points,edges,centers):
    
    x = points[:,0].flatten()
    y = points[:,1].flatten()

    plt.figure(figsize=(10,8))
    plt.plot(centers[:,0],centers[:,1],'o',color ='g')
    plt.plot(x[edges.T], y[edges.T], linestyle='-', color='y',markerfacecolor='red', marker='o') 
 