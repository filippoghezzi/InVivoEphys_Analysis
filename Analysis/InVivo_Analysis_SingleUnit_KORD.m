close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData.mat')
data.PSTHvisualOpto=[];
data.PSTHlaser=[];
data.responseVisualOpto=[];
data.responseLaser=[];
data1=data;
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData_Part2.mat')
data=[data1;data];
folderFigures='C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\Conferences\11th Annual Oxford Neuroscience Symposium\MyInVivo';
