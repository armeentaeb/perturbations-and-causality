close all
clear all

%% setting up global parameters
global Sigma; % consists of covariances across environments
Sigma = cell(7,1);
global learn_rate
global nvec
learn_rate = 0.001;
global Bstar;
global latent
latent = 1;
global perturbed_latent; perturbed_latent = 1;
global psi_max; psi_max = 1;
global rank_est; rank_est = 2;
global R;
global ratio_BIC; ratio_BIC = 1;


%% setting up working directories
temp = pwd;
root = temp(1:36);

addpath(strcat(root,'/cvx'));
cvx_startup
cvx_setup
addpath(strcat(root,'/Solvers'));


Bstar = load(strcat(root,'/files/population_Bstar'));


%% Adjust this for the particular 
link = strcat(root,'/files/Fig2_zeta5_h1_psi1_latent/');

num_obs = 300;
nvec = [num_obs 5 5 5 5 5 5];
p = 10;

R = ones(p,length(nvec));

% size of the intervental data
iterv = [64,16,4,2,1];
FP_final = zeros(p^2,length(iterv)); TP_final = zeros(p^2,length(iterv));
for j =1:length(iterv)
    
    link_in = strcat(link,string(iterv(j)*5),'_samples/');
    nvec = [num_obs 5*iterv(j)*ones(1,6)];
    
    
    FP_tot = zeros(p^2,1);
    TP_tot = zeros(p^2,1);
    FP_tot_singleparam = zeros(p^2,1);
    TP_tot_singleparam = zeros(p^2,1);
    
    % averaging across 10 independent trials
  
    
    for iter = 1:10
        link_final = strcat(link_in,'equivalent_DAGS_GES',string(iter));
        
        % obtaining all the data
        tot_data = load(link_final);
        
        % DAGS from pooled GES
        equivalent_DAGs = cell(size(tot_data,1)/p,1);
        for l = 1:size(tot_data,1)/p
            equivalent_DAGs{l} = tot_data(1+(l-1)*p:p*l,:);
        end
        link_final = strcat(link_in,'overall_data',string(iter));
        tot_data = load(link_final);
        data_obs = tot_data(1:floor(num_obs*8/10),:);
       
        % the covariance matrices across environments
        Sigma{1} = cov(data_obs);
        for iter_2 = 1:length(Sigma)-1
            Sigma{iter_2+1} = cov(tot_data(num_obs+5*iterv(j)*(iter_2-1)+1:num_obs+5*iterv(j)*(iter_2-1)+5*iterv(j),:));
        end
        
        % scoring DAG
        [optimal_Bv1,Gamma,d0,scorev1] = score_all_dags_set(equivalent_DAGs);
        
        
        
        %% performing backward deletion
        [val1,ind1] = min(scorev1);
        optimal_B_est_orig = optimal_Bv1{ind1};
        
        
        
        M = optimal_Bv1{ind1};
        [row col] = find(abs(optimal_Bv1{ind1})>0.001);
        vec = zeros(length(row),1);
        for i = 1:length(row)
            vec(i) = M(row(i),col(i));
        end
        [val ind] = sort(abs(vec),'ascend');
        row = row(ind);
        col = col(ind);
        
        scorev2 = 0;
        score1 = scorev1(ind1);
        score_val = zeros(length(row),1);
        score_val(1) = scorev1(ind1);
        result = M;
        DAG{1} =  M;
        
        
        cov_test = cov(tot_data(1+floor(num_obs*8/10):num_obs,:));
       
        clear test
        precision= (eye(p)-optimal_Bv1{ind1})'*(diag(d0{1})+Gamma{1}*Gamma{1}')^(-1)*(eye(p)-optimal_Bv1{ind1});
        test(1) = -log(det(precision))+trace(precision*cov_test);
        
        for l = 1:length(row)
            clear equivalent_DAGs2;
            equivalent_DAGs2 = [];
            M(row(l),col(l))= 0;
            equivalent_DAGs2{1} = M;
            
            [optimal_Bv2,Gamma,d0,scorev2] = score_all_dags_set(equivalent_DAGs2');
            M = optimal_Bv2{1};
            scorev2
            score_val(l+1) = scorev2;
            DAG{l+1} =  M;
            if scorev2> score1
                %break;
            else
                result = M;
                score1 = scorev2;
            end
             precision= (eye(p)-optimal_Bv2{1})'*(diag(d0{1})+Gamma{1}*Gamma{1}')^(-1)*(eye(p)-optimal_Bv2{1});
             test(l+1) = -log(det(precision))+trace(precision*cov_test);
        end
        
        % choosing the DAG with best test error
        [val ind] = min(test);
        result = DAG{ind};
        B_opt = result;
        B_opt = abs(B_opt)>10^(-3);
        % finding the error in connectivity matrices
        [FP,TP] = evaluate_error_connectivity(B_opt);
        FP_tot = FP_tot+FP; TP_tot = TP_tot+TP;
    end
    FP_tot = FP_tot/10;
    TP_tot = TP_tot/10;
    FP_final(:,j) = FP_tot; TP_final(:,j) = TP_tot;
end