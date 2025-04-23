addpath(genpath(pwd))
warning('off')

%% problem and field distribution for academic test example 
savetoMat = false;
deg = 2;
subs = 10;
f = 150;
%academicTestExample;
clear all

%% problem and field distributions for planar coil example
tic;
deg = 2;
subs = 10;
f = 150;
%planarCoilExample;
toc
clear all

%% problem and field distributions for inductively coupled example
tic;
deg = 2;
subs = 10;
f = 150;
%inductiveCouplingExample;
toc
clear all

%% problem and field distributions for planar coil using the unstabilized original formulation
tic;
deg = 2;
subs = 10;
f = 150;
%noGauge_planarCoil;
toc
clear all

%% Generation of condition number over frequency plots
% problem setup has to be changed in the method! 
% currently: academic
deg = 2;
subs = 2;
freqs = logspace(-6, 14, 11);
%conditionOverFreq;
clear all

%% convergence study for quantitative validation using original formulation
% currently still has to be run multiple times for different sigma values
tic;
convergenceStudyOriginal;
toc
clear all

%% convergence study for quantitative validation using stabilized formulation
% currently still has to be run multiple times for different sigma values
tic;
%convergenceStudyStabilization;
toc
clear all





