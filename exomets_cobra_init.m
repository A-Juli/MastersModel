% Automated start of day process for setting up COBRA and loading the base
% version of the model.

initCobraToolbox

% modelDirectory='C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO' ;
% modelFileName_base = 'iCHOv1_S_final.xml';
modelDirectory='C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO' ;
modelFileName_base = 'GSM_mouse_trimmed_BDMFA.xml';


modelFileName_base = [modelDirectory filesep modelFileName_base];
% model_iCHOv1_S_Base = readCbModel(modelFileName_base);
model_BDMFA_Base = readCbModel(modelFileName_base);
% compSymbolList=model_iCHOv1_S_Base.comps;
% compNameList=model_iCHOv1_S_Base.compNames;
compSymbolList=model_BDMFA_Base.comps;
compNameList=model_BDMFA_Base.compNames;

% 
% fprintf('\n\n');
% fprintf('Model Parameters Cheatsheet \n');
% fprintf('Units are mmol per gDW per hr \n');
% fprintf('Compartment IDs: \n');
% fprintf('c  = Cytosol \n');
% fprintf('r  = Endoplasmic Reticulum \n');
% fprintf('e  = Extracellular Space \n');
% fprintf('g  = Golgi Apparatus \n');
% fprintf('x  = Peroxisome \n');
% fprintf('im = Intermembrane Space of Mitochondria \n');
% fprintf('m  = Mitochondria \n');
% fprintf('l  = Lysosome \n');
% fprintf('n  = Nucleus \n');

