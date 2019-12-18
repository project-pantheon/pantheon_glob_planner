##!/usr/bin/env python3

import os
import numpy as np

def getHomeDir():

    return os.path.expanduser("~")

def getCurrWDir():

    return os.getcwd()


def writeFile(file_name, object_to_write, specifics):

    file_path = getCurrWDir() + "/miscellaneous/"

    f = open(file_path + file_name,"w")

    f.write( str(object_to_write) )
    
    f.close()

#    np.savetxt(file_path + file_name, object_to_write, fmt = specifics, delimiter=',')

def generateEucTSPFile(namefile,num_tree,points):
    
    file_path = getCurrWDir() + "/inputs/"
    file_name = "EUC_Grid_th_" + str(inBetweenRows) + ".tsp"
    
    f = open(file_path + file_name,"w")
    f.write( "NAME: %s" % file_name )
    f.write( "\nTYPE: TSP" )
    f.write( "\nCOMMENT: waypoints" )
    f.write( "\nDIMENSION: %d" % n_points )
    f.write( "\nEDGE_WEIGHT_TYPE: EUC_2D" )
    f.write( "\nNODE_COORD_SECTION\n" )
    
    
    for i in range(points):
        s = str(i+1) + " "+ str(points[i][0])+" "+ str(points[i][1])+ "\n"
        f.write(s)
        
    f.write("EOF")
    f.close()

def generateTSPFileFromAdj(inBetweenRows, n_points, W, AdjList):

    file_path = getCurrWDir() + "/inputs/"
    file_name = "Grid" + str(inBetweenRows) + ".tsp"

    f = open(file_path + file_name,"w")

    f.write( "NAME: %s" % file_name )
    f.write( "\nTYPE: TSP" )
    f.write( "\nCOMMENT: waypoints" )
    f.write( "\nDIMENSION: %d" % n_points )
    f.write( "\nEDGE_WEIGHT_TYPE: EXPLICIT" )
    f.write( "\nEDGE_WEIGHT_FORMAT: FULL_MATRIX" )

    #not taken into account ---> CHECK ISSUE
    f.write( "\nEDGE_DATA_FORMAT: ADJ_LIST" )
    
    f.write( "\nEDGE_DATA_SECTION" )
    for i,j in AdjList.items():
        oneString = str(j).replace(',', '')
        s = "\n" + str(i) + " "+ oneString.strip('[]') + " -1 "
        f.write(s)

    f.write("-1")
    f.write( "\nEDGE_WEIGHT_SECTION\n" )
    
    
    for i in range(n_points):
        for j in range(n_points):
            s =  str(W[i][j]) + " "
            f.write(s)
        s = "\n"
        f.write(s)

 
    f.write("EOF")
    f.close()

    print("TSP file correctly written")


def writeOutCome(one_stop_time,inBetweenRows, n_rows, n_cols, opt_dist , solution, comp_time):

    file_path = getCurrWDir() + "/results/"
    file_name = "Sol_th_" + str(inBetweenRows) + "_" + str(n_rows) + "x" + str(n_cols) + ".txt"

    f = open( file_path + file_name, "w")    
    
    f.write( "TSP Results for a %s x %s grid considering a inBetweenRows %s" % (str(n_rows), str(n_cols), str(inBetweenRows) ) )
    f.write( "\nOptimal distance (Euclidian sense): %s [m]" % opt_dist )
    f.write( "\nTSP total computational time: %s [sec]" %comp_time )
    f.write( "\nTotal number of stops: %s" % str(solution.shape[0]) )
    f.write( "\nTotal stop time results: %s [sec]" % str(solution.shape[0]*one_stop_time) )
    f.write( "\nThe optimal nodes sequence results:\n" )
    f.write( str(solution) )
    f.close()

    print("Results correctly written")
