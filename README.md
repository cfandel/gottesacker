# Groundwater modeling for the Gottesacker karst aquifer
*Chloé Fandel, 2019*
<br>*University of Arizona Dept. of Hydrology & Atmospheric Sciences*
<br>*Karlsruhe Institute of Technology Institut für Angewandte Geowissenschaften*

This project uses Python to model groundwater flow in karst, by linking a 3D geologic model (using GemPy), a conduit network model (using SKS), and a pipe flow model (using SWMM).
Many possible conduit networks are generated in SKS based on the geologic model, then flow is routed through each proposed conduit network using SWMM and the predicted spring discharge time series are compared to real data. 
<br>Project details: https://www.agw.kit.edu/mitarbeiter_10064.php

<br>*GemPy:* open-source, Python-based 3-D structural geological modeling software from CGRE Aachen. https://github.com/cgre-aachen/gempy
<br>*SKS:* Stochastic Karst Simulator - pseudo-genetic 3D stochastic modeling of karst aquifers from University of Neuchatel. https://doc.rero.ch/record/31792
<br>*SWMM:* Storm Water Management Model - runoff pipe flow model from US EPA. https://www.epa.gov/water-research/storm-water-management-model-swmm

### Primary project notebooks:

- ensemble_generator.ipynb: this is the bulk of the project code to generate an ensemble of models
- ensemble_viewer.ipynb: this is the plotting & analysis once the ensemble has already been generated

### Modules being worked on:

*For GemPy:*
- mapping.py: additional visualization functions to crop the 3D model based on the land surface or irregularly-shaped boundaries, to export and import gslib files, and to plot custom diagonal cross-sections
- manager.py: creates correctly formatted GemPy input files by selecting subsets of data from a dataframe with all available data in it

*For SKS:*
- sksin.py: create input files for SKS based on a template, and run SKS as a subprocess from Python. 

*For SWMM:*
- swmmpy.py: create input files for SWMM based on a template, and run SWMM as a subprocess from Python. 


 
