%% This is a demo for the paper "A Statistical Graphical Model of the California Reservoir Network
%% Written by Armeen Taeb, California Institute of Technology, December 2016
close all
clear all
clc

%%  Define globa variables for solvers
global GM; % Graphical modeling
global LVGM; % Latent variable graphical modeling
global Conditional_LVGM; % Conditional latent variable graphical modeling
GM = 1;
LVGM = 2;
Conditional_LVGM = 3;
global covariatesName; % cell containing name of the covariates
global resInfo
global TrainTest
global root
global Sigma
global nvec
global learn_rate
global perturbed_latent
global latent
global psi_max
global rank_est
global R



% get appropriate directory
temp = pwd;
root = temp(1:end);
root1 = temp(1:end-22);


% addpath to the code
addpath(strcat(root1,'/cvx'));
cvx_startup
cvx_setup

addpath(strcat(root1,'/Solvers/'));


addpath(root1)
addpath(strcat(root,'/reservoir_data'));



% read reservoir data
[DataY] = xlsread('ReservoirsData');
DataY = DataY(:,2:end);
selRes = [1:9,11:13,15:16,18:29,31:49,51:60]; % These are reservoirs without huge chunk of missing information
DataY = DataY(:,selRes);

[resInfoM,resInfoT] = xlsread('reservoir_summary');

% numeric fields
resInfo.lat = resInfoM(selRes,1);
resInfo.lon = resInfoM(selRes,2);
resInfo.elev = resInfoM(selRes,3);
resInfo.cap = resInfoM(selRes,4);
resInfo.drain = resInfoM(selRes,5);
resInfo.zoneind = resInfoM(selRes,7);
resInfo.basinind = resInfoM(selRes,9);
resInfo.streamind = resInfoM(selRes,11);


% text fields
resInfo.name = resInfoT(selRes+1,1);
resInfo.dam = resInfoT(selRes+1,2);
resInfo.lake = resInfoT(selRes+1,3);
resInfo.hydro = resInfoT(selRes+1,4);
resInfo.basinname = resInfoT(selRes+1,5);
resInfo.zonename = resInfoT(selRes+1,6);
resInfo.basinname = resInfoT(selRes+1,7);
resInfo.streamname = resInfoT(selRes+1,8);


p = 10;
numYears = floor(size(DataY,1)/365);
numMonths = numYears*12 + floor((size(DataY,1)/365-numYears)*12);
% Obtain monthly average reservoir volumes
DataYM = AverageMonthly(DataY,numMonths);
[val,ind] = sort(resInfo.cap,'descend');
DataYM = DataYM(:,ind(1:p));

timePeriod = 1:size(DataYM,1); % start in January 2004 since snowpack data starts



%% Create training and validation data
%validation observations from January 2003 - December 2003, January
%2013 - November 2015
% Training and Validation Data
TestInd = [];
TrainInd = setdiff(1:size(DataYM,1),TestInd);
TestY = DataYM(TestInd,:); % reservoirs validation data
TrainY =  DataYM(TrainInd,:); % reservoirs training data
TrainX = DataYM(TrainInd,:); % covariates training data
TestX = DataYM(TestInd,:); % covariates validation data


[TrainYa, TestYa, avgMonthY,varY] = preprocessingData(TrainY,TestY,TrainInd,TestInd);
data = TrainYa;
tot_env = 4;
%env_ind{1}=[1:48,85:10*12]; %normal
env_ind{1}=[1:48,85:8*12];
env_ind{2} = [48+1:12*5,10*12+1:12*11]; %abnormally dry: 2007, 2013
env_ind{3} = [12*5+1:7*12]; %moderought drought: 2008-2009
env_ind{4} = [12*11+1:12*13];

Test_Data = data(8*12+1:10*12,:);

A = load('water_GES2');
equivalent_DAGs  = cell(length(A)/p,1);
for i = 1:length(A)/p
    equivalent_DAGs{i} =  A((i-1)*p+1:(i)*p,:);
end



Sigma_train = cell(tot_env,1);
for i = 1:tot_env
    Sigma_train{i} =  cov(data(env_ind{i},:));
    nvec_train(i) = length(env_ind{i});
end

Sigma = Sigma_train;
nvec = nvec_train;
R = ones(p,length(nvec_train));


learn_rate = 0.00001;
perturbed_latent = 1;
latent = 1;
psi_max = 2;
rank_est = 2;


B_opt_wl = learning_and_prediction(equivalent_DAGs,Test_Data);



% without latent variables
learn_rate = 0.00001;
perturbed_latent = 1;
latent = 0;
psi_max = 2;
rank_est = 2;


B_opt_wol = learning_and_prediction(equivalent_DAGs,Test_Data);


