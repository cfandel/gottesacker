%@** SKSmain/options_karst
% PURPOSE
% This file is used to set the options for the main script "sks.m"
% Placeholder strings for certain variables can be replaced using the python wrapper module sksin.py 
% All commas have been replaced with semicolons and all blank lins are preceded by a % 
%  in this version to enable python compatibility.
% See comments below to understand the options.
% 
% VERSION LOG
% v01: - wild version unusable by normal human beings
% v02: - new version usable by normal human beings
% v03: - modified by Louis Fournier
% v04: - new options
%      - clean version
%      - translated in english
% v05: - new features such as:
%           * system Id (there might be several karst systems within the same model) 
%           * different fractures families depending on geological formation
% 
% AUTHOR
% Andrea Borghi 
% modified by Chloé Fandel
% 
% CREATION DATE
% June 2009
% last update: August 2016
% python template: 2019
%
% DISCLAIMER
% This software is an experimental research prototype. The software is
% provided by the copyright holders and contributors "as is" without
% warranty of any kind; express or implied; including; but not limited to;
% the implied warranties of merchantability and fitness for a particular
% purpose and non infringement. In no event shall the authors or copyright
% owner or contributors be liable for any claim; any direct; indirect;
% incidental; special; exemplary; or consequential damages (including; but
% not limited to; procurement of substitute goods or services; loss of use;
% data; or profits; or business interruption) however caused and on any
% theory of liability; whether in contract; strict liability; or tort
% (including negligence or otherwise) arising in any way out of the use of
% or in connection with this software; or other dealing in the software
% even if advised of the possibility of such damage.
%@
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL VARIABLES DECLARATION (DO NOT EDIT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FMvalFaults;
global FMvalFrac;
global multiplicatorConduits;
global codeAquiclude;
global codeAquifere;
global codeOut;
global kindOfParticleTrack; 
global reRunFMA_YoN;
global verbosity;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDIT OPTIONS IN THE FOLLOWING SECTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL OUTPUT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% level of "verbosity" of sks. I.e. the quantity of outputs in the command
% line terminal during execution.
% from 0 (no output) to 3 (full output) - suggested value: 1 (2 and 3 are more 
% intended for debug reasons).
verbosity=1;
%
% simulation basename: all the output files will start with this string
nomSimul='output_SKS'; 
%
% create a separate output directory (1) or not (0)
outputDirYoN=1; 
%
% name of the separate output directory
outputDir=cd;
%
% output MAT-file with the full workspace
saveWorkspaceYoN=0;
%
% output VTK files for visualization
VTK_YoN=1;
%
% output GSLIB files
outputGslib=0;
%
% output ASCII files for mesh and conduits
saveAsciiMeshFiles=1; % - 0: no output
                      % - 1 ascii files only for mesh 1D
                      % - 2 ascii files only for mesh 3D
                      % - 3 ascii files for both 1D and 3D meshes
%
% output ASCII LIN and DAT files for Feflow
saveAsciFeflowFiles_YoN=1; % 0 -> no ; 1 -> yes; set to yes if simplifying network afterwards
%
% Save option files (usefull to remember what has been done)
logOptions_YoN=0;
%
% output a VTK file of the model DTM (with model grid resolution)
DTM=0;  
%
% save the inlets of the model
saveStrPts_YoN=0; 
%
% save a gslib file with the fractures
saveFracturesYoN=0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MESH DIMENSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% origin of the grid
x0=578287.5;    %xmin from DEM raster
y0=5240062.5;   %ymin from DEM raster
z0=800.0;       %zmin from DEM raster
%
% number of elements in each direction
ni=200; % in y dimension: i.e. nb of rows (yres from DEM)
nj=251; % in x dimension: i.e. nb of columns (xres from DEM)
nk=52;  % in z dimension: i.e. nb of layers (z res from GemPy)
%
% size of elements
dx=50;      %pixel width
dy=50;      %pixel height
dz=27.5;      %dz = (zmax-zmin)/z res
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FAST MARCHING SPEEDS FOR THE DIFFERENT MODEL FACIES (very important)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remember: the higher the value  the higher is the influence of this
% formation. 
%
codeOut=0.0011;                 % outside the topography FMA value (should be close to 0)
codeAquifere=2;            % aquifer FMA value
codeAquiclude=1;           % aquiclude FMA value (should be smaller than codeAquifere)
%
FMvalFrac = 20;             % fractures maximal FMA value
FMvalFaults = 1000;           % faults FMA value
FMvalHorizons = 6;         % inception horizons FMA value
%
multiplicatorConduits = 2;  % multiplicator for the conduits: FMvalConduits = FMvalFrac * multiplicatorConduits
%should be set so that conduits are more conductive than the aquifer
%
useCumulativeFractures_YoN=0;   % if 0: milieu(fractures>=1)=FMvalFrac    
                                % if 1: milieu(fractures>=1)=fractures(fractures>= 1)*(FMvalFrac-codeAquifere)/max(max(max(fractures)))+codeAquifere;
                                % this mean that if useCumulativeFractures_YoN is set  the FMA
                                % value for the fractures is taken between the value of the aquifer 
                                % and the value "FMValFrac"
%                         
% if reRunFMA_YoN is set: the fast marching will be re-runned if the walker
% tries to go to a cell that is above its starting altitude. This adds
% consistency to the generated conduits.
reRunFMA_YoN = 1; % 0 :no - usually for debug reasons
                  % 1: yes - common use
%                  
% Simultate first the vadose conduits toward a base level (the altitude of
% the spring) and secondly toward the spring.
step1 = 1; % 0: no -  usually for debug reasons
           % 1: yes -  common use
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INLETS OF THE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% type of inlets. Random -  load from file etc.
startPtsType = 2; % - 0 = random inlet points. The user must define xMinStart 
                          %   yMinStart xMaxStart  and yMaxStart as bounding box
                          %   and also the "nb_RandomInlets" variable (see below).
                          % - 1 = load from an xy file -  uses "importanceFactor" to 
                          %   define the proportions of inlets in each iteration
                          %   (see below)
                          % - 2 = random inlets + load from file. Uses nb_RandomInlets
                          %   to know the number of random inlets and then uses
                          %   importanceFactor as in case 1.
                          % - 3 = load from file (xyc) and class the inlets
                          %   by catchment surface (3th column of the file)
                          % - 4 = mix of 1 and 3 -> random inlets + loaded inlets
                          %   and class them by catchment. Again uses nb_RandomInlets to
                          %   known the number of random inlets and
                          %   importanceFactor for the hierarchy.
                          % - 5 = classified by system ID (input = [x y z id])
%                 
% startingPtsInFile : the xy file containing the coordinates of the inlets
startingPtsInFile='input/inlets.txt';
%
% For Gottesacker: use known inlets from Zhao's model + random inlets
%
% randomly permute the inlet coordinates (add variability to the models)
permuteCoord=1;
%
% "surrounding box" for the random inlets (lower left and upper right
% corners)
xMinStart=583667.252;   %lower left corner coordinates
xMaxStart=586108.038;   %upper right corner coordinates
yMinStart=5245286.976;
yMaxStart=5247000.156;
%
%xMinStart=x0;          
%xMaxStart=xMinStart*nj*dx;   %extend to entire model
%yMinStart=y0;
%yMaxStart=yMinStart*ni*dy;
%
% number of random inlets. 
nb_RandomInlets=0; 
%
importanceFactor = [1 1 1 1 1 ];   % importanceFactor. It is a vector:
                            % - its length correspond to the number of iterations
                            % - the number that are stored correspond to the relative
                            %   number of loaded inlets that will be in the several iterations.
                            %   this number is defined as the current element over the sum
                            %   of importanceFactor
                            %
                            % EXAMPLE:
                            % if importanceFactor = [1 3 6]   -> sum(importanceFactor)=10
                            % -> there will be 3 iterations
                            % -> the first will contain the first 1/10 of the total inlets
                            % -> the second will contain the next 3/10 of the total inlets
                            % -> and the third the remaining inlets
%                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTLETS OF THE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% springInputXYZFile: the coordinates input file for the karst springs
springInputXYZFile='input/springs.txt'; 
%
% order of computation
permuteSprings=1;        % randomly permute spring order (useless if less than 2 springs)
alternateSpringsYoN = 1; % 0 : all the springs are considered in only one FMA iteration
                         % 1 : every spring is iterated sequentially  i.e
                         % multi-connections of springs and sinkholes
%
% karstSystemConnectivity -> {i} = spring class "i" 
%the array in "i" element correspond to the inlet classes that are connected to these springs|, 
%if {} is empty then the default behavior is used
karstSystemConnectivity = {};
%
% karstSystemSimulationOrder -> which sub system to simulate and in which
% order. This allow to deactivate useless systems
% karstSystemSimulationOrder = [1 2 3 8]; % SW systems
% karstSystemSimulationOrder = 0; -> deactivate this feature
karstSystemSimulationOrder = 0; %[1];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEOLOGICAL MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% type of geological model: 'normal'  'geomodeller' or 'petrel'
% In the case of a geomodeller or petrel output file  has to be reshaped in
% a particular way. 
% It typeOfGeology is "uniform"  it will be created by SKS as a homogeneous
% field (filled with codeAquifer)
% NOTE: what is 'normal'?
typeOfGeology='geomodeller';  %use geomodeller formate for GemPy files
%
% input file for geology     %if in different folder|, use '../foldername/filename.gslib'             
geoleFileName = 'input/geo.gslib'; %can have either only lith or lith and orientation
%
% geologyIdArray: array with the geological indexes in the "geoleFileName"
% that need to be considered.
%i.e. the formation numbers from GemPy: no faults
    %NaN             0  ##this is for inactive cells (outside model bounds)
    %Garschella      1                 
    %Schrattenkalk   2  ##this is the aquifer               
    %Drusberg        3                 
    %basement        4 
%
% with faults:
    %%NaN            0  ##this is for inactive cells (outside model bounds)
    %faultB1         1
    %faultA2         2
    %Garschella      3                 
    %Schrattenkalk   4  ##this is the aquifer               
    %Drusberg        5                 
    %basement        6 
%
geologyIdArray=[1 2 3 4];
%
% geologyFMArray: the corresponding FMA values that have to be given to the
% formations contained in "geologyIdArray". Both arrays must have the same
% size -  if some formations are not listed in these arrays but are in the
% geological model  a value of "codeOut" will be assigned to them.
% assign 'codeAquiclude' to aquicludes  and 'codeAquifere' to aquifer (aka
% karstifiable formation)
geologyFMArray=[codeAquiclude codeAquifere codeAquiclude  codeAquiclude];
%
inception_YoN=0; % 1 = random inception horizons 
                 % 2 = all the inception horizons
nb_inception=1;  % if inception_YoN==1 must provide the number of randomly-selected inception horizons
inceptionIdArray=[2]; % similar to "geologyIdArray" -  tells which ids in the geological model can be considered as inception horizons
%
% cut the model borders using a polygon
cutModelByPolygon= 0; 
% CSV file containing the polygon for the model limit (x y z - no header)
cutPolygonCsvFile='./inputFiles/limitPolygon.csv'; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FRACTURES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fracturesYoN = 2;      % 0 : no fractures
                          % 1 : load fractures from gslib file
                          % 2 : use the built-in fracture generator
%
% input gslib file for fractures (if fracturesYoN == 1)
fracturesFileName='input_sks/faults_for_SKS.gslib';
%
% wich geological formations need fractures
geologyFractureIDs=[2];
%
% total number of fractures to be generated (for each geology)
chooseNFrac=1500;  
%
% strike: code 99999 for random
% Orientations BAGET : N0; N45 et N115 (+90° dans l'algorithme)
% Gottesacker: N310  N20 (add 90 for SKS directions)
%so 310 + 90 = 40& 20 + 90 = 110
%
strikeMinFamiliesCell={[35 105 ]};
strikeMaxFamiliesCell={[45 115 ]};
%
% dip for all the families: code 99999 for random
dipMinFamiliesCell={[80 80 ]}; %[0 80 80 80 80 80 99999]; 
dipMaxFamiliesCell={[90 90 ]};
%
% similar to "importanceFactor": proportion of fractures in each family.
% again it has to be defined for each geological formation
importanceFactorFamFracCell={[8 2 ]};  
%
% type of fracture family: extensive (1)  compressive (-1) or neutral (0)
typeOfFractureRegimeCell={[0 0 0 0 0 0 0 ]};
%
% cell array containing the vectors with the lenghts for fractures
% (one vector for each fracture family) - > look into adding
lll=[25:0.1:500];
fractureLengthCell={lll};   
%
% type of fracture length distribution ('UFM' is only option currently available)
fractureLengthDistributionType={'UFM'};
%
% to force fracture to be conditional to inlet points
fracturesOnInletPoints_YoN=0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISCRETE FAULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use faults or not 
faultsYoN = 0; % 0 : no faults
               % 1 : load discrete faults from gslib file
               % 2 : load discrete faults from xyz file (petrel format)
               % 3 : load discrete faults from simple xyz files
%              
% input gslib file for faults (binary matrix  0=void  1=faults)
faultsFileName='input_SKS/faults_for_SKS.gslib'; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPLEX CONDUIT PATHS (plugs in the medium to force diverging conduits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connectDeconnect_YoN = 0; % 0 no plugs:
                          % 1 "balls" plugs  "rBalls" has to be defined and
                          % 2 gaussian plugs  paramsGauss has to be defined
nPlugsMax=0;    % maximal number of plugs between one iteration and the other
                  % NB: the true value is randomly taken between nPlugsMax/2 and
                  % nPlugsMax.
%
% min and max radius for "balls" plugs
rBalls=[0.01 10];
%
% paramsGauss = [lx ly lz  mu et sigma2  quantile]
% where : lx ly and lz are the correlation length in x y and z dimensions
% mu is the mean of the distribution
% sigma2 is the variance
% quantile is the quantile threshold of cut
paramsGauss=[300 300 300 0 1 .95];  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLOW SIMULATION FILES AND MESHING (unstable -  may need code modifications)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MODFLOW Finite-Differences simulator
modflowYoN=0; % create files for MODFLOW (1) or not (0)
%
% GroundWater -  Finite-Element simulator 
% if onlyCutTheMesh is set -  the FE mesh is cut outside "cutPolygon"
onlyCutTheMesh = 0;
%
use1Dconduits = 0; % mesh 1D pipes in the finite element mesh (very time consuming)
% NB: to create valid input files for GW  at least one of the 2 above
% options must be set to 1.
%
GW_YoN=0; % create files for GW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL  DO NOT EDIT THE FOLLOWING SECTION (but leave it...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% save also vtk files for each generated fracture (long and heavy)
FaireVTK_Fractures=0; 
%
only1DunstrtuctMesh=0; % 0: mesh3D avec / sans conduits suivant use1Dmesh
                       % 1: mesh elements finis seulement 1D
%
% cut the model underneath the deeper aquifer (reduce CPU times for flow
% simulation by deactivating useless aquitard elements)
cutBelow=0;
%
% used in antamina branch: see reclassStuff.m
nbConduitsOrdersForOutput=7;
%
% niveaux de base automatiques
baseLevelYoN = 0; % si on veut que l'algo cherche des sources par rapport a un niveau de base  utile pour aquiferes cotiers
baseLevel = [3700 4000]; % vecteur contenant dans l'ordre les niveaux de base
nbSprings = 1; % nombre de sources qu'on veut par niveau de base -  attention c'est defini automatiquement pour le cas ou les sources sont connues -  ctd baseLevelYoN=0
%
options_karst_file=which('options_karst');
options_gw_file=which('options_create_GW_stuff');
%
FMorGW = 1; % switch between FM (1) or steady-state GW simulation (2) as dist_map
            % GW: still to do for 3D and step1 !!!! do not use!!!!
%
kindOfParticleTrack=2; % 1 = connectivite a 9 noeuds: 2 = connectivite a 6 noeuds
%
densityFromOrientation=0; % 0 or 1: if yes (1) it will put more fractures in
                          % the points of the model with biggest intensity of deformation
%
classifyPropFieldYoN = 1; % 0 : propField will be output with real values
                          % 1 : propField will be output classified
%  
typeDeSimul='humanVersion';
%
% epikarst layer (not really relevant)
epikarstZone=0; % create an epikarst layer on the top of model (1) or not (0)
nbrLayersEpi=1; % number of elements of the epikarst layer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF EXPERIMENTAL SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
