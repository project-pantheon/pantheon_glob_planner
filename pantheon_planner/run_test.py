#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from concorde.tsp import TSPSolver
import functions as fc
import plots_functions as pf
import file_writer as fw
import time
import sys
from shapely.geometry import LineString
from shapely.geometry import Point
from termcolor import colored
from collections import defaultdict
from concorde.tests.data_utils import get_dataset_path

np.set_printoptions(threshold=sys.maxsize)


def computeTsp(n_trees,dist,points):   
    
    generateTSPFile("nodes.tsp",n_trees,dist)
    solver = TSPSolver.from_tspfile("/home/majo/libs/pyconcorde/pantheon_planner/inputs/nodes.tsp")

    start_time = time.time()
    solution = solver.solve()
    print( colored("--- %s seconds ---" % (time.time() - start_time)) )
    print( colored("Outcome TSP:","green") )
    print(solution.found_tour)
    print( colored("With optimal value:","green") )
    print(solution.optimal_value/100)
 
    sol = solution.tour

    print(sol)
#    edges = fc.formatEdges(sol)
#    fo = computeFo(dist,edges) #don't need to calculate the distance since is stored in the optimal_value
#    print(fo)
#    return fo,sol
    return solution.optimal_value,sol


def computeAdjTsp(filename):   
 
    solver = TSPSolver.from_tspfile("/home/majo/libs/pyconcorde/pantheon_planner/inputs/" + filename)

    start_time = time.time()
    solution = solver.solve()
    eval_time = time.time() - start_time
    
    print( colored("--- %s seconds ---" % eval_time) )

    print( colored("Outcome TSP:","red") )
    print(solution.found_tour)
    print( colored("With optimal value [m]:","red") )
    print(solution.optimal_value/100)
 
    sol = solution.tour
    print(solution)
    print(sol)
#    edges = fc.formatEdges(sol)
#    fo = computeFo(dist,edges) #don't need to calculate the distance since is stored in the optimal_value
#    print(fo)
#    return fo,sol
    return solution.optimal_value,sol,eval_time

def distMatrix(n_trees, points,centers): #ASYMMETRIC MATRIX
    dist = np.zeros(shape=[n_trees*4,n_trees*4],dtype=int)    
    for i in range(n_trees*4):
        for j in range(n_trees*4):
           if ((i%4 == 0) and (j == i+3)) or ((i%4 == 1)and(j ==i+1)):
                dist[i][j] = 1000000
                dist[j][i] = 1000000
           else: 
               dist[i][j] = np.round(np.sqrt(np.power(points[i][0]-points[j][0],2)+np.power(points[i][1]-points[j][1],2)))
                
    return dist    

def computeAdjMatrix(n_rows,n_cols,inBetweenRows):

    #given the n° of trees we have distribution based on inBetweenRows 
    if ( inBetweenRows == 1):
        MAP_DIM = ( n_rows + 1 ) * ( n_cols + 1 )
        BLOCK_SIZE = ( n_cols + 1 );
        N_BLOCKS = ( n_rows + 1) ;
        v1 = np.ones( (MAP_DIM - BLOCK_SIZE) );

        COLS_connections = np.add( np.diag(v1,BLOCK_SIZE), np.diag(v1,-BLOCK_SIZE) );

        v2 = np.ones( (BLOCK_SIZE-1) );

        BLOCK = np.add( np.diag(v2,1), np.diag(v2,-1) );

        ROWS_connections = np.zeros((MAP_DIM,MAP_DIM))

        for i in range(N_BLOCKS):
            ROWS_connections[ i*BLOCK_SIZE : (i+1) * BLOCK_SIZE, i*BLOCK_SIZE : (i+1) * BLOCK_SIZE] = BLOCK;
    
        Adj_matrix = np.add( COLS_connections, ROWS_connections );
    
    else:
        MAP_DIM = ( ( n_rows + 1) * n_cols ) + ( ( n_cols + 1) * n_rows )
        
        v1 = np.ones( (MAP_DIM - n_cols) );
        v2 = np.ones( (MAP_DIM - (n_cols + 1) ) );


        for i in range(1,n_rows+1):
            v1[ n_cols + (i-1)*(2*n_cols +1) ] = 0;

        for i in range(1,n_rows):
            v2[2*n_cols+ (i-1)*(2*n_cols +1)] = 0

        upper = np.add( np.diag(v1,n_cols) , np.diag(v2, n_cols + 1) ) 

        lower = np.transpose(upper) 

        Adj_matrix = np.add(upper,lower)


    return Adj_matrix


def fromMatrix2List(A):

    graph = defaultdict(list)
    edges = set()

#    graph2 = defaultdict(list)
#    edges2 = set()

    for i, v in enumerate(A, 1):
        for j, u in enumerate(v, 1):
            if u != 0 and frozenset([i, j]) not in edges:
                edges.add(frozenset([i, j]))
                graph[i].append(j)   # put graph[i].append({j: u}) if you need also to have the weight 
				     # in this case the TSP file writing has to be modified
#    for i, v in enumerate(A, 1):
#        for j, u in enumerate(v, 1):
#            if u != 0 and frozenset([i-1, j-1]) not in edges2:
#                edges2.add(frozenset([i-1, j-1]))
#                graph2[i-1].append(j-1)   # put graph[i].append({j: u}) if you need also to have the weight 
				     # in this case the TSP file writing has to be modified

    return graph


def computeWeightMatrix(Adj,points):

    dist = np.zeros(shape=[Adj.shape[0],Adj.shape[0]],dtype=int)    
    for i in range(Adj.shape[0]):
        for j in range(Adj.shape[0]):
                dist[i][j] = np.round(np.sqrt(np.power(points[i][0]-points[j][0],2)+np.power(points[i][1]-points[j][1],2)))


    #include the tree constraints directly by weight matrix

    return dist


def addConstraints(W, n_rows, n_cols, inBetweenRows):

    if ( inBetweenRows ): 
        for i in range(1,W.shape[0]-(n_cols)):
            if ( i%(n_cols+1) == 1 ):
                W[i-1][i+n_cols+1] = 999
                W[i+n_cols+1][i-1] = 999
            elif ( i%(n_cols+1) == 0 ):
                W[i-1][i+n_cols-1] = 999
                W[i+n_cols-1][i-1] = 999
            else:
                W[i-1][i+n_cols-1] = 999
                W[i-1][i+n_cols+1] = 999
                W[i+n_cols-1][i-1] = 999
                W[i+n_cols+1][i-1] = 999
    else:
        for i in range(n_rows):
            for j in range(n_cols):
                W[i*(2*n_cols+1) + j][(i+1) * (2*n_cols+1) +j] = 999
                W[(i+1) * (2*n_cols+1) +j][i*(2*n_cols+1) + j] = 999

        for i in range(n_cols):
            for j in range(n_rows):
                W[n_cols + j*(2*n_cols+1) + i ][n_cols +j*(2*n_cols+1) +i+1] = 999
                W[n_cols +j*(2*n_cols+1) +i+1][n_cols + j*(2*n_cols+1) + i ] = 999


    return W


def generateTSPFile(namefile,n_tree,distMatrix):
    
    f = open("inputs/"+namefile,"w")
    f.write( "NAME: %s" % namefile )
    f.write( "\nTYPE: TSP" )
    f.write( "\nCOMMENT: waypoints" )
    f.write( "\nDIMENSION: %d" % (n_tree*4) )
    f.write( "\nEDGE_WEIGHT_TYPE: EXPLICIT" )
    f.write( "\nEDGE_WEIGHT_FORMAT: FULL_MATRIX" )
    f.write( "\nEDGE_WEIGHT_SECTION\n" )
    
    
    for i in range(n_tree*4):
        for j in range(n_tree*4):
            s =  str(distMatrix[i][j])+" "
            f.write(s)
        s = "\n"
        f.write(s)
    f.write("EOF")
    f.close()



def main(R=5,n_rows=4, n_cols=2): #REMARK: n_rows*n_cols MUST BE ODD, OTHERWISE HAMILTONIAN PATH DOESN'T EXIST


#PARAMETERS
#R = 1000 radious of trees    
    inBetweenRows = 0 # where we want to take points around trees. REMARK: with n_rows = n_cols = EVEN no HAMILTONIAN path can be found.
                      #for small rows and cols it finds still a (wrong) solution, while for large rows and cols it fails 
    obstacles_enabled = 1 #trees are considered in the weighted graph
    molt = 1000
    tree_radious = 0.3
    dist_row = 5.0#*molt
    dist_col = 5.0#*molt
    n_trees = n_cols*n_rows

    est_single_stop_time = 3*60 #sec


#    My solution
   
    trees_pos = fc.generateTrees(n_rows,n_cols,dist_row,dist_col)
    
#    pf.plotMapTrees(trees_pos,R)
    if (inBetweenRows == 0):

        capture_pts = fc.generateCapturePts(trees_pos,n_rows,n_cols, dist_row,dist_col) #ok

    else:
        capture_pts = fc.generateinBetweenCapturePts(n_rows,n_cols, dist_row,dist_col) #ok
#    pf.plotCapturePts(capture_pts)

    print(capture_pts.size)

 

    pf.plotTreesAndCaptures(trees_pos,capture_pts,tree_radious)

    Adj_matrix = computeAdjMatrix(n_rows,n_cols,inBetweenRows)
    
    W = computeWeightMatrix(Adj_matrix,capture_pts)

    if ( obstacles_enabled ):
        W = addConstraints(W,n_rows,n_cols,inBetweenRows)

    Adj_list = fromMatrix2List(Adj_matrix)
    
    fw.generateTSPFileFromAdj(inBetweenRows,capture_pts.shape[0], W*100, Adj_list) #TSP solver takes only int numbers 

    #TSP solver assumption: complete graph.
    opt_dist , soluti, comp_time = computeAdjTsp("Grid" + str(inBetweenRows) + ".tsp")

    edges = fc.formatEdges(soluti)
    
    pf.plotSolutionGraph(capture_pts,edges,trees_pos,"SolGrid" + str(inBetweenRows) + str(n_rows) + "x" + str(n_cols) + ".txt")

    fw.writeOutCome(est_single_stop_time,inBetweenRows, n_rows, n_cols, opt_dist/100 , soluti, comp_time)

    sys.exit()

#CHECK IF AN HAMILTONIAN PATH IS POSSIBLE 

#FILE FORMAT:

#NAME : alb1000
#COMMENT : Hamiltonian cycle problem (Erbacci) 
#TYPE : HCP
#DIMENSION : 1000
#EDGE_DATA_FORMAT : EDGE_LIST
#EDGE_DATA_SECTION
#  1000    593
#  1000    456


   ##### ALESSANDRO ######
    centers = fc.generateCenters(n_rows,n_cols,dist_row,dist_col)

    points_ref = fc.generateRefFrame(centers,R)

    pf.plotCapturePts(points_ref)

    if (inBetweenRows == 0): #in-between trees
        thetas = np.zeros(shape=(1,n_trees))#np.pi/4#+(np.pi/2)*np.random.rand(1,n_trees)
        points = fc.computePoints(thetas, points_ref, centers)
        dist = distMatrix(n_trees,points,centers)

        fo,sol = computeTsp(n_trees,dist,points)
        edges = fc.formatEdges(sol)
        pf.plotGraph(points,edges,centers,n_cols,n_rows)
        pf.plotPoints(points,edges)
        print(colored("inBetweenRows = 0, distance [km] : " +str(fo/molt),"cyan"))
    else:

        thetas45 = np.zeros(shape=(1,n_trees))+3*np.pi/4#+(np.pi/2)*np.random.rand(1,n_trees)
        points45 = fc.computePoints(thetas45, points_ref, centers)
        dist45 = distMatrix(n_trees,points45,centers)
        fo45,sol45 = computeTsp(n_trees,dist45,points45)
        edges45 = fc.formatEdges(sol45)
    
        pf.plotGraph(points45,edges45,centers,n_cols,n_rows)
        pf.plotPoints(points45,edges45)
        print(colored("inBetweenRows = 45 distance [km]: " +str(fo45/molt),"green"))
#    time.sleep(0.5)
  

#    pf.plotWithoutLines(points,centers,n_cols,n_rows)
#    
   
    sys.exit()


    
if __name__ == '__main__':  
    main()