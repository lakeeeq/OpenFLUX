clearvars
clc
%%%run this script to create files produced previously by JAVA
%%%files are created at the same path as the input model
addpath OFfunctions %can be added to MATLAB startup

OF = OpenFLUX;
OF.modelFileName = 'OF_2009_projectFolder/modelInput_sample.txt';%metabolic model
OF.isDynamic = false;
OF.version = '2009';

OF.buildModel;