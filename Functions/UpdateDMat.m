%Align EEG times independently of DMAI program for updating DMAI files
%This program will overwrite original file, make a backup!
%Load Dmat File
clear; 
[file, pf] = uigetfile('.mat','Select a DMAT file'); 
[EEG, pe] = uigetfile('.csv','Select EEG File'); 
load([pf file]); 
[~,~,raw] = xlsread([pe EEG]); 
X = raw(2:end,1:20); 
Y = cell2mat(X); 

%Replace DMAT EEG with selected EEG
D{1}.eegData = Y; 
D{1}.eegFilename = EEG; 

%Set current directory to that of the file to be replaced
cd(pf); 
%Overwrite EEG in D matrix of selected DMAT file
%save(file,'D','FILES','param1','param2'); 
save(file,'D', 'param1', 'param2', 'FILES'); 

