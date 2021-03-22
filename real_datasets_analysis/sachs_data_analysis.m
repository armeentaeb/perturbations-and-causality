close all
clear all
clc



% get appropriate directory
temp = pwd;
root = temp(1:end);
root1 = temp(1:end-22);

% addpath to the code
addpath(root1)
addpath(strcat(root,'/cvx'));
cvx_startup
cvx_setup


% names of the proteins
Names = {'b2camp.xls', 'cd3cd28+aktinhib.xls','cd3cd28+g0076.xls',...
     'cd3cd28+ly.xls','cd3cd28+psitect.xls','cd3cd28+u0126.xls',...
     'cd3cd28.xls','pma.xls'};

global Sigma;
global nvec;
p = 11;
Select_env = setdiff(1:8,[7]);

Sigma = cell(length(Select_env)+1,1);
data = log(xlsread(strcat(root,'/Sachs_data/',Names{7})));
Sigma_train{1} = cov(data(1:floor(size(data,1)*8/10),:));
nvec_train(1) = floor(size(xlsread(strcat(root,'/Sachs_data/',Names{7})),1)*8/10);
Sigma_test{1} = cov(data(1+floor(size(data,1)*9/10):end,:));
Test_Data = data(1+floor(size(data,1)*9/10):end,:);
nvec_test(1) = size(data,1)-nvec_train(1);


l = 2;
train_env = setdiff(Select_env,[]);
for i = 1:length(train_env )
    data = log(xlsread(strcat(root,'/Sachs_data/',Names{i})));
    Sigma_train{l} = cov(data(1:floor(size(data,1)*10/10),:));
    nvec_train(l) = size(data,1);
    l = l+1;
end

Sigma = Sigma_train;
nvec = nvec_train;
global R;
R = ones(p,length(nvec_train));

%% now generate equivalent DAGs based on the truth
equivalent_DAGs = cell(336,1);
for i = 1:336
       equivalent_DAGs{i} = load(strcat(root,'/Sachs_data/Equivalent_DAGS/equivalent_DAGS_',string(i)));
end

global learn_rate
global perturbed_latent
global latent
global psi_max
global rank_est
learn_rate = 0.0005;
perturbed_latent = 1;
latent = 1;
psi_max = 0.5;
rank_est = 2;

[B_opt_wl] = learning_and_prediction(equivalent_DAGs,Test_Data);



% without latent variables
learn_rate = 0.0005;
perturbed_latent = 1;
latent = 0;
psi_max = 0.5;
rank_est = 2;

[B_opt_wo] = learning_and_prediction(equivalent_DAGs,Test_Data);




















