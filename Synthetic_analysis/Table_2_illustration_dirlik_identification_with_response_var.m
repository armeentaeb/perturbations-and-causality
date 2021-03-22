close all
clear all

global Sigma; % consists of covariances across environments
global env_num;
Sigma = cell(5,1);
global learn_rate
global nvec
learn_rate = 0.01;
global Bstar;
global latent
latent = 1;
global perturbed_latent; perturbed_latent = 1;
global psi_max; psi_max = 2;
global rank_est; rank_est = 2;
global R;
global ratio_BIC; ratio_BIC = 1;




temp = pwd;
root = temp(1:36);

addpath(strcat(root,'/cvx'));
cvx_startup
cvx_setup
addpath(strcat(root,'/Solvers'));

Bstar = load(strcat(root,'/files/population_Bstar'));
Bstar(10,3) = -0.7;
Bstar(10,4) = -0.7;
Bstar(:,10) = 0;
Bstar(7:9,10) = -0.7;

%% Adjust this for the particular 



num_obs = 1000;
nvec = [num_obs 5 5 5 5];
p = 10;
R = ones(p,length(nvec));


iterv = [60,200];
obs = 5*iterv;
num_env = 5;

psi_max = 2;
%[FD_strong_latent_intervention,TD_strong_latent_intervention] = DirectLikelihood(link,p,num_env,iterv,obs);

link =  strcat(root,'/files/comparison_no_intervention_DAG/');
[FD_no_intervention_DAG,TD_no_intervention_DAG] = DirectLikelihood(link,p,num_env,iterv,obs);

link =  strcat(root,'/files/comparison_target_intervention_DAG/');
[FD_target_intervention_DAG, TD_target_intervention_DAG] = DirectLikelihood(link,p,num_env,iterv,obs);

link =  strcat(root,'/files/comparison_latent_intervention_DAG/');
[FD_latent_intervention_DAG,TD_latent_intervention_DAG] = DirectLikelihood(link,p,num_env,iterv,obs);

link =  strcat(root,'/files/comparison_all_intervention_DAG/');
[FD_all_intervention_DAG,TD_all_intervention_DAG] = DirectLikelihood(link,p,num_env,iterv,obs);


