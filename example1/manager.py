'''Python tools to manage GemPy-SKS-SWMM interactions & files.
Chloé Fandel, 2018.
cfandel@email.arizona.edu

Miscellaneous process management tools.

See accompanying notebook for examples of how to use. 
'''

#################################################################################################
import numpy as np
import pandas as pd
import os 
from copy import copy

#################################################################################################

def create_gp_input(allpts, unitnames, orientations_filename, interfaces_filename):
        
        '''Creates correctly formatted GemPy input files by selecting data from a dataframe with all possible points in it.
        Dataframe must have a column 'use' with values 'y'/'n' indicating whether to use each data point or not.
        Inputs:
        allpts:                pandas df of all possible data pts to choose from. 
                               Must have columns: sourceID, ID, type, X, Y, Z, formation, azimuth, dip, polarity, 
                               xsec surface, xsec bottom, thickness, zsurf (DEM), surface elev diff (DEM - Goldscheider), 
                               use, interface, overturned, Source, notes
        unitnames:             list of strings indicating formation names to be used in this model
        orientations_filename: string or path indicating name of orientations file to be created.
        interfaces_filename:   string or path indicating name of interfaces file to be created.
        
        Outputs:
        orientations.csv:      csv file of selected orientation points in correct format for GemPy (columns: X,Y,Z,azimuth,dip,polarity,formation)
        interfaces.csv:        csv file of selected interface points in correct format for GemPy   (columns: X,Y,Z,formation)
        '''
        
        allpts = allpts[allpts.use=='y']                      #select only data points flagged as yes to use (can change flag based on whatever criteria desired)
        allpts = allpts[allpts.formation.isin(unitnames)]     #select only data points for units in the list of unit names being used (defined above)
        orientations = allpts[allpts.type=='orientation']     #select orientation pts (need 1 pt per fm)
        interfaces   = allpts[(allpts.type=='interface') | (allpts.type=='xsec')]     #select interface pts (need 2 pts per fm)
        
        #check number of occurrences of each formation (need at least 1 pts per fm for orientations, and 2 pts per fm for interfaces)
        if np.any(orientations.formation.value_counts().values<1) or len(np.unique(orientations.formation.values))!=len(unitnames):              
            print('STOP! Need at least 1 point per formation in orientations data. At least one formation is missing data:\n', orientations.formation.value_counts())
        if np.any(interfaces.formation.value_counts().values<2) or len(np.unique(interfaces.formation.values))!=len(unitnames):               
            print('STOP! Need at least 2 points per formation in interfaces data. At least one formation is missing data:\n', interfaces.formation.value_counts())

        orientations.sort_values(['order','X','Y'], inplace=True)                #sort by formation then location
        try:
            orientations.drop(columns=['ID','type','order','use','fault','overturned','source','erosion'], inplace=True)  #drop unneeded columns
        except:
            orientations.drop(columns=['ID','type','order','use'], inplace=True)
        colnames = ['X', 'Y', 'Z', 'azimuth', 'dip', 'polarity', 'formation']    #list of column names being kept, in the correct order
        orientations = orientations.reindex(columns=colnames)                    #reorder the column names
        
        interfaces.sort_values(['order','formation','X'], inplace=True)          #sort by formation
        try:
            interfaces.drop(columns=['ID','type','order','use','fault','overturned','source','erosion'], inplace=True)       #drop unneeded columns
        except:
            interfaces.drop(columns=['ID','type','order','use'], inplace=True)
        colnames = ['X', 'Y', 'Z', 'formation']                                  #list of column names being kept, in the correct order
        interfaces = interfaces.reindex(columns=colnames)                        #reorder the column names
        
        orientations.to_csv(orientations_filename, index=False)                  #export to csv
        interfaces.to_csv(  interfaces_filename,   index=False)                  #export to csv
