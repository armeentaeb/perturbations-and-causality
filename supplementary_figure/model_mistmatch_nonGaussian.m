close all
clear all

%% setting up global parameters
global Sigma; % consists of covariances across environments
Sigma = cell(7,1);
global learn_rate
global nvec
learn_rate = 0.01;
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
link = strcat(root,'/files/mismatch_nonGaussian/');

num_obs = 300;
p = 10;
nvec = [1000 5 5 5 5 5 5];
R = ones(p,length(nvec));


equivalent_DAGs = cell(9*3,1);
for i = 1:9
    equivalent_DAGs{3*i-2} = load(strcat('equivalent_DAGS_',string(i)));
    B = equivalent_DAGs{3*i-2};
    ind = find(B == 1);
    isDAG = 0;
    while ~isDAG
        Bnew = B;
        t = randperm(length(ind));
        Bnew(ind(t(1:2))) = 0;
        sumv = 0;
        for j = 1:p
            sumv = sumv + trace(Bnew^j);
        end
        isDAG = (sumv == 0);
    end
    equivalent_DAGs{3*i-1} = Bnew;
   
    
    ind = find(B == 0);
    isDAG = 0;
    while ~isDAG
        Bnew = B;
        t = randperm(length(ind));
        Bnew(ind(t(1:2))) = 1;
        sumv = 0;
        for j = 1:p
            sumv = sumv + trace(Bnew^j);
        end
        isDAG = (sumv == 0);
    end
    equivalent_DAGs{3*i} = Bnew;
end

iterv = [64,16,4];
FP_final = zeros(p^2,length(iterv)); TP_final = zeros(p^2,length(iterv)); 
for j =1:length(iterv)
    link_in = strcat(link,string(iterv(j)*5),'_samples/');
    nvec = [num_obs 5*iterv(j)*ones(1,6)];

    FP_tot = zeros(p^2,1);
    TP_tot = zeros(p^2,1);
    FP_tot_singleparam = zeros(p^2,1);
    TP_tot_singleparam = zeros(p^2,1);

    for iter = 1:10
         link_final = strcat(link_in,'equivalent_DAGS_GES',string(iter));
        tot_data = load(link_final);
        equivalent_DAGs = cell(size(tot_data,1)/p,1);
        for l = 1:size(tot_data,1)/p
            equivalent_DAGs{l} = tot_data(1+(l-1)*p:p*l,:);
        end
        link_final = strcat(link_in,'overall_data',string(iter));
        tot_data = load(link_final);
        data_obs = tot_data(1:num_obs,:);
        Sigma{1} = cov(data_obs(1:floor(size(data_obs)*8/10),:));
        cov_test = cov(tot_data(1+floor(num_obs*8/10):num_obs,:));

        for iter_2 = 1:length(Sigma)-1
            Sigma{iter_2+1} = cov(tot_data(num_obs+5*iterv(j)*(iter_2-1)+1:num_obs+5*iterv(j)*(iter_2-1)+5*iterv(j),:));
        end
        
        [optimal_Bv1,Gamma,d0,scorev1] = score_all_dags_set(equivalent_DAGs);
        
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
        
        
       clear test
       precision= (eye(p)-optimal_Bv1{ind1})'*(diag(d0{1}))^(-1)*(eye(p)-optimal_Bv1{ind1});
       test(1) = -log(det(precision))+trace(precision*cov_test);
        
        scorev2 = 0;
        score1 = scorev1(ind1);
        result = M;
        for l = 1:length(row)
            clear equivalent_DAGs2;
            equivalent_DAGs2 = [];
            M(row(l),col(l))= 0;
            equivalent_DAGs2{1} = M;
            
            [optimal_Bv2,Gamma,d0,scorev2] = score_all_dags_set(equivalent_DAGs2');
            M = optimal_Bv2{1};
            scorev2
            score_val(l) = scorev2;
            DAG{l} =  M;
            if scorev2> score1
                %break;
            else
                result = M;
                score1 = scorev2;
            end
            
             precision= (eye(p)-optimal_Bv2{1})'*(diag(d0{1})+Gamma{1}*Gamma{1}')^(-1)*(eye(p)-optimal_Bv2{1});
             test(l+1) = -log(det(precision))+trace(precision*cov_test);
        end
        
        [val ind] = min(test);
        B_opt = DAG{ind};
        [FP,TP] = evaluate_error_connectivity(B_opt);
        FP_tot = FP_tot+FP; TP_tot = TP_tot+TP;
               
    end    
    FP_tot = FP_tot/10;
    TP_tot = TP_tot/10;
    FP_final(:,j) = FP_tot; TP_final(:,j) = TP_tot;  
end
