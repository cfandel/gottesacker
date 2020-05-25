'''Python tools to interact with SKS (Stochastic Karst Simulator).
Chloé Fandel, 2019.
cfandel@email.arizona.edu

Uses a template options_karst.m input file to enable editing parameters in different sections of the input file. 
Can also run SKS as a subprocess from Python.

See accompanying notebook for examples of how to use. 
See SKS documentation for descriptions of parameters.

Notes: 
- Not all parameters are enabled for editing yet, and many defaults are not yet modifiable, 
but the basic structure of the code should make it easy to add functions in the same pattern.
- the template file CANNOT have any commas (otherwise the csv writer thinks they are escape characters)
- the template file CANNOT have any blank lines - use a % symbol at the beginning of each blank line (otherwise the csv writer has problems)
- placeholders in template file must be unique
- order of placeholder list and data list must be the same
- common error: the template filename isn't correct, so the insert_data() function can't find the section where it's supposed to insert the data

SKS created by Andrea Borghi, 2009 (last update 2016)
SKS Copyright: Sclumberger Water Services, Paris, & University of Neuchatel
SKS Disclaimer: 
% This software is an experimental research prototype. The software is
% provided by the copyright holders and contributors "as is" without
% warranty of any kind, express or implied, including, but not limited to,
% the implied warranties of merchantability and fitness for a particular
% purpose and non infringement. In no event shall the authors or copyright
% owner or contributors be liable for any claim, any direct, indirect,
% incidental, special, exemplary, or consequential damages (including, but
% not limited to, procurement of substitute goods or services; loss of use,
% data, or profits; or business interruption) however caused and on any
% theory of liability, whether in contract, strict liability, or tort
% (including negligence or otherwise) arising in any way out of the use of
% or in connection with this software, or other dealing in the software
% even if advised of the possibility of such damage.
'''


########################################################################################
import subprocess as sb
import pandas as pd
import matplotlib.pyplot as plt
import swmmpy as sp


#########################################################################################

def import_template(template_filename='options_karst_template.m'):
    
    '''Imports the SKS options_karst_template.m file into a pandas dataframe.
    Inputs:
    template_filename: filename/path for the options_karst_template.m file
                       this must be a text file correctly formatted for SKS, with placeholder strings in the sections to be edited.
                       To create a sample template file, use the options_karst.m file from the SKS documentation, and replace the data
                       in the MESH DIMENSIONS/origin of the grid section with the strings 'xmin','ymin','zmin'
    Outputs:
    template:          pandas dataframe of the template file, with one column, and a line for each line of text.'''
    
    template = pd.read_csv(template_filename, header=None, skip_blank_lines=False, delimiter='\n', encoding='ISO-8859-1')  #import template file to dataframe
    
    return template


#########################################################################################

def insert_data(template, placeholders, data, show=False):
    '''Inserts data from a pandas dataframe into a template dataframe by replacing the section placeholder in the template.
    The template file and data file must be read in first. Best to make a copy of the template.
    
    Inputs:
    template:     pandas dataframe with one column, and one line per line of text in the template input file.
                  NOTE: each section that will be edited MUST have a section placeholder string at the desired insert location
    placeholder:  list of strings corresponding to the placeholder string to be replaced
    data:         list of parameter values to insert (must be same length and order as placeholder list)
    show:         True/False. If True, prints the updated section
    
    Outputs:
    new:          pandas dataframe with one column, and one line per line of text in the desired output file, 
                  with the new data inserted where the placeholder was in the original file
    '''
    
    for i in range(len(placeholders)):     #loop over placeholders
        placeholder = placeholders[i]      #select the current placeholder from list
        value       = data[i]              #select the value to insert
        if type(value) == type([]):        #if the value is a list
            value = ''.join(str(x)+' ' for x in value)
            value = '[' + value + ']'
        else:
            value = str(value)             #convert to string
        #print('replacing ', placeholder, 'with ', value)
        mask = template[0].str.contains(placeholder)   #get a True/False array for which rows contain the placeholder
        mask[mask.isna()] = False                      #replace nans with False
        line = template[mask][0]                       #select the line with the placeholder
        #print('line', line.index.values[0], ': \t', line.values[0])
        new = line.str.replace(placeholder, value)  #replace string in current line with the desired value
        #print('newline: \t', new.values[0])     
        template[mask] = new                           #replace original line in template with updated line
        if show==True:
            print('updated line: \t', template[mask][0].values[0], '\n') #print updated line if requested
            
    return template


#########################################################################################

def write_input(placeholders, data, template_filename='options_karst_template.m', inputfile='options_karst.m'):
    
    '''Write a SWMM input file using a template.inp file with placeholders in each sections,
    and pandas dataframes holding the data to be inserted for each section. 
    The placeholders will be replaced with data.
    Section order must be the same in placeholders and data.
    
    Inputs:
    inputfile:           string for name of input file to write (defaults to options_karst.m). 
    placeholders:        list of placeholder strings to be replaced with data. 
    data:                list of parameters to be inserted 
    template_filename:   string for name of template file (defaults to 'options_karst_template.m') - MUST NOT HAVE COMMAS!!!!
    
    Outputs:
    options_karst.m:    SKS input file
    
    Note: the template options_karst.m file cannot include any commas (because these are treated as escape characters and cause problems)
    '''
    
    template = import_template(template_filename)                       #import template file from txt
    template = insert_data(template, placeholders, data)                #replace placeholder string with data
    template.to_csv(inputfile, header=False, index=False, quoting=3, escapechar='|')    #write dataframe to .inp text file with specified name
    #template.to_csv(inputfile, header=False, index=False, quoting=3)    #write dataframe to .inp text file with specified name

    
#########################################################################################    
    
def run_SKS(sks_path, matlab='C:/Program Files/MATLAB/R2017b/bin/matlab.exe'):
    
    '''Run SKS from MATLAB as a subprocess. All SKS settings are controlled by the options_karst.m file and the run_SKS.m file.
    options_karst.m can be modified using the write_input() file and a template file.
    See SKS documentation for explanation of variables.
    Note: Paths to SKS folder and to matlab exe must be strings with forward slashes and no leading r'''
    
    #SKS options:
    sks       = 'runSKS'                           #name of script to run
    options   = '-nosplash -nodesktop -wait'        #optional: set run options  (nosplash? nodesktop means MATLAB won't open a new desktop window, wait means Python will wait until MATLAB is done beore continuing (needs to be paired with p.wait() after sb.Popen))
    has_args  = False                               #set whether the MATLAB script needs arguments (i.e. is it a function?)

    #Set function string:
    #Structure:  """path_to_exe optional_arguments -r "cd(fullfile('path_to_folder')), script_name, exit" """
    #Example:    """C:/Program Files/MATLAB/R2017b/bin/matlab.exe -r "cd(fullfile('C:/Users/Chloe/Desktop/PhD/SKS/')), run_SKS, exit" """
    fun =  """{} {} -r "cd(fullfile('{}')), {}, exit" """.format(matlab, options, sks_path, sks)  #create function string that tells subprocess what to do

    #Run MATLAB:
    print('running SKS in MATLAB...')
    p = sb.Popen(fun)                     #open the subprocess & run the MATLAB script 
    p.wait()                              #wait until MATLAB is done before proceeding (this needs to be paired with -wait in options?)
    print('done')                         #if the run is successful, an output file names a.mat should appear in the folder with the MATLAB scripts


#########################################################################################    

def plot_network(nodes_file, links_file, spring_file, elev_file=None, dim=3, labels=False):
    
    '''Imports simplified network info (nodes and links) from SKS outputs to pandas dataframes for nodes and links,
    then plots the network map. 
    
    Inputs:
    nodes_file:  text file with two (or three) columns, for X, Y (and Z) coordinates of simplified node locations
    links_file:  text file with two columns, for upstream node and downstream node of each link
    spring_file: text file with three columns, for X,Y,Z coordinates of springs (aka outlets)
    elev_file:   text file with one column, for elevation of each node (only if 2D and elev added later)
    dim:         number of dimensions (2 or 3)
    labels:      True/False, whether to display node labels
    
    Outputs:
    nodes:       pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    links:       pandas df with all links, and columns: InNode, OutNode, Name
    '''
    
    if dim == 2:
        nodes         = pd.read_csv(nodes_file, header=None, names=['X','Y'], delim_whitespace=True) #import node data from SKS txt output file
        nodes['Z']    = pd.read_csv(elev_file, header=None, names=['Z'])                    #import node elev data from SKS and add it as a new column
    
    if dim == 3:
        nodes         = pd.read_csv(nodes_file, header=None, names=['X','Y','Z'], delim_whitespace=True) #import node data from SKS txt output file
        
    nodes['Name'] = nodes.index                                                         #add a column for node name (based on index)
    springs       = pd.read_csv(spring_file, delim_whitespace=True)                     #import spring location data
    
    links         = pd.read_csv(links_file, header=None, names=['InNode','OutNode'], delim_whitespace=True) #import link data from SKS txt output file (simplified)
    links         = links - 1                                                           #convert node indices to Python 0-based  indices
    links['InNode']  = [int(node) for node in links.InNode]                             #convert to integer
    links['OutNode'] = [int(node) for node in links.OutNode]                            #convert to integer
    links['dz']      = [nodes.loc[links.loc[i].InNode].Z - nodes.loc[links.loc[i].OutNode].Z for i in links.index]  #add new column for elevation change between two nodes
    links = links[links.dz.values>0]                                                    #select only conduits with a downhill slope
    links.drop(labels=['dz'], axis='columns', inplace=True)                             #drop  unneeded columns
    links.reset_index(drop=True,inplace=True)                                           #reset indexing to start from zero
    links['Name'] = links.index                                                         #add a column for link name (use row index)
        
    fromX = nodes.X.loc[links.InNode]                                                   #calculate coordinates for link start and end points
    fromY = nodes.Y.loc[links.InNode]
    toX   = nodes.X.loc[links.OutNode]
    toY   = nodes.Y.loc[links.OutNode]
    
    if dim==3:
        fromZ = nodes.Z.loc[links.InNode]
        toZ   = nodes.Z.loc[links.OutNode]
    
    if dim==2:
        f,ax = plt.subplots()                                                                #create figure & axis objects
        ax.scatter(nodes.X, nodes.Y, color='k', s=8)                                         #plot nodes
        for ind in links.index:                                                              #loop over link indices
            plt.plot((fromX.iloc[ind],toX.iloc[ind]),(fromY.iloc[ind],toY.iloc[ind]), c='k') #plot links
        if labels==True:
            for ind in nodes.index:                                                          #loop over node indices
                plt.annotate(ind,xy=(nodes.X[ind],nodes.Y[ind]+100))                         #label nodes with index
    if dim==3:
        f = plt.figure(figsize=(10,10))                                                 #create empty figure
        ax = f.add_subplot(111, projection='3d')                                        #add 3D subplot axes (requires Axes3D module from mpl_toolkits.mplot3d)
        ax.scatter(nodes.X, nodes.Y, nodes.Z, color='k', s=8)                           #plot nodes
        for ind in links.index:                                                         #loop over link indices
            ax.plot((fromX.iloc[ind],toX.iloc[ind]),(fromY.iloc[ind],toY.iloc[ind]),(fromZ.iloc[ind],toZ.iloc[ind]), c='k')    #plot links
        if labels==True:
            for ind in nodes.index:                                                          #loop over node indices
                ax.text(nodes.X[ind],nodes.Y[ind],nodes.Z[ind],  '%s' % (str(ind)))          #label nodes with index
                
    ax.set_aspect('equal')
    plt.title('SKS nodes & junctions')
    
    