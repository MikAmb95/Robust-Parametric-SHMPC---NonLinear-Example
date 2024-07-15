
clc,close all,clear all
% Get the current directory
currentFolder = pwd;

% Add the subfolder and all its subdirectories to the MATLAB path
addpath(genpath(fullfile(currentFolder, 'YalmipFolder')));
addpath(genpath(fullfile(currentFolder, 'MPT3')));
addpath(genpath(fullfile(currentFolder, 'TUMcps-CORA-537aec3')));

addpath(genpath(fullfile(currentFolder, 'Results')));
addpath(genpath(fullfile(currentFolder, 'ControlExtraFunction')));
addpath(genpath(fullfile(currentFolder, 'helperFunctions')));

disp('First Script: pRSHMPC L=2')
pause(1)
run('SpaceCraftAttitudeDyn_p_RHMPC_L2.m')
clc
disp('Second Script: pRSHMPC L=10')
pause(1)
run('SpaceCraftAttitudeDyn_p_RHMPC_L10.m')
clc
disp('Third Script: pRSHMPC L=20')
pause(1)
run('SpaceCraftAttitudeDyn_p_RHMPC_L20.m')
clc
disp('Foruth Script: RSHMPC')
pause(1)
run('SpaceCraftAttitudeDyn_RHMPC.m')

clc
disp('--END--- press enter to print results')
pause()
run('printResults.m')
