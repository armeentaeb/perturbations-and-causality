close all
clear all
global Sigma; % consists of covariances across environments
Sigma = cell(7,1);
global learn_rate
global nvec
learn_rate = 0.01;
global Bstar;
global latent
latent = 1;
global perturbed_latent; perturbed_latent = 1;
global psi_max; psi_max = 0.5;
global rank_est; rank_est = 3;
global R;
global ratio_BIC
ratio_BIC = 1;



temp = pwd;
root = temp(1:36);

addpath(strcat(root,'/cvx'));
cvx_startup
cvx_setup

Bstar = load(strcat(root,'/files/population_Bstar'));
addpath(strcat(root,'/Solvers'));


num_obs = 1000;
nvec = [num_obs 5 5 5 5 5 5];
p = 10;
R = ones(p,length(nvec));

link = strcat(root,'/files/');
link_in = link;
for l = 1:8
    link_final = strcat(link_in,'equivalent_DAGS_',string(l));
    M = load(link_final);
    ind = find(M == 0);
    isDAG = 0;
    while ~isDAG
        M = load(link_final);
        M = M+eye(p);
        ind = find(M == 0);
        t = randperm(length(ind));
        M(ind(t(1:5))) = 1;
        M = M-eye(p);
        sumv = 0;
        for j = 1:p
            sumv = sumv + trace(M^j);
        end
        isDAG = (sumv == 0);
    end
    equivalent_DAGs{l} =M;
    
end

M = Bstar; M(5,10) = 1;M(8,10) = 1;M(5,3)=1;
equivalent_DAGs{9} = M;

for l = 9:9
    
    
    isDAG = 0;
    while ~isDAG
        M = equivalent_DAGs{l};
        M = M+eye(p);
        ind = find(M == 0);
        t = randperm(length(ind));
        M(ind(t(1:2))) = 1;
        M = M-eye(p);
        sumv = 0;
        for j = 1:p
            sumv = sumv + trace(M^j);
        end
        isDAG = (sumv == 0);
    end
    %M = load(link_final);
    equivalent_DAGs{l} =M;
    
end

equivalent_DAGs = equivalent_DAGs';


link =  strcat(root,'/files/comparison_to_two_stage_non_dense/');
iterv = [200];
Sigma_deconf = cell(5,10,2);
Sigma_normal = cell(5,10,2);

for j =1:length(iterv)
    link_in = strcat(link,string(iterv(j)*5),'_samples/');
    for iter = 1:10
        
        
        link_final = strcat(link_in,'overall_data',string(iter));
        tot_data = load(link_final);
        data_obs = tot_data(1:num_obs,:);
        
        
        [precision,graph_model] = deconfound(data_obs);
                
        Sigma_normal{1,iter,j} = cov(data_obs(1:floor(size(data_obs,1)*8/10),:));
        Sigma_deconf{1,iter,j} = graph_model^(-1);
        for iter_2 = 1:4
            [precision,graph_model] = deconfound(tot_data(num_obs+5*iterv(j)*(iter_2-1)+1:num_obs+5*iterv(j)*(iter_2-1)+5*iterv(j),:));
            Sigma_deconf{iter_2+1,iter,j} = graph_model^(-1);
            Sigma_normal{iter_2+1,iter,j} = cov(tot_data(num_obs+5*iterv(j)*(iter_2-1)+1:num_obs+5*iterv(j)*(iter_2-1)+5*iterv(j),:));
        end
    end
end








FP_final = zeros(p^2,length(iterv)); TP_final = zeros(p^2,length(iterv));
FP_final_deconfound = zeros(p^2,length(iterv)); TP_final_deconfound = zeros(p^2,length(iterv));
FP_final_backshift = zeros(p^2,length(iterv)); TP_final_backshift = zeros(p^2,length(iterv));

FP_final_thr = zeros(p^2,length(iterv)); TP_final_thr = zeros(p^2,length(iterv));
FP_final_deconfound_thr = zeros(p^2,length(iterv)); TP_final_deconfound_thr = zeros(p^2,length(iterv));

for j =1:length(iterv)
    link_in = strcat(link,string(iterv(j)*5),'_samples/');
    nvec = [num_obs 5*iterv(j)*ones(1,4)];
    
    
    FP_tot = zeros(p^2,1);
    TP_tot = zeros(p^2,1);
    FP_tot_deconfound = zeros(p^2,1);
    TP_tot_deconfound = zeros(p^2,1);
    
    FP_tot_backshift  = zeros(p^2,1);
    TP_tot_backshift  = zeros(p^2,1);
    
    
    FP_tot_thr = zeros(p^2,1);
    TP_tot_thr = zeros(p^2,1);
    FP_tot_deconfound_thr = zeros(p^2,1);
    TP_tot_deconfound_thr = zeros(p^2,1);
    
    for iter = 1:10
        
        
        link_final = strcat(link_in,'overall_data',string(iter));
        tot_data = load(link_final);
        data_obs = tot_data(1:num_obs,:);
        
        
        
        Sigma = Sigma_normal(:,iter,j);
        latent = 1;
        perturbed_latent = 1;
        [optimal_Bv1,Gamma,d0,scorev1] = score_all_dags_set(equivalent_DAGs);
        [val1,ind1] = min(scorev1);
        cov_test = cov(tot_data(1+floor(num_obs*8/10):num_obs,:));
        
        
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
                [val ind] = min(score_val);
                M = DAG{ind};
     
        B_opt = M;
        B_opt = abs(B_opt)>10^(-3);
        [FP,TP] = evaluate_error_connectivity(B_opt);
        
        
        
        Sigma = Sigma_deconf(:,iter,j);
        latent = 0;
        perturbed_latent = 0;
        
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
        
        scorev2 = 0;
        score1 = scorev1(ind1);
        score_val = zeros(length(row),1);
        score_val(1) = scorev1(ind1);
        result = M;
        DAG{1} =  M;
        
                clear test
                precision= (eye(p)-optimal_Bv1{ind1})'*(diag(d0{1}))^(-1)*(eye(p)-optimal_Bv1{ind1});
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
                [val ind] = min(score_val);
        
               M = DAG{ind};
  
        B_opt = M;
        B_opt = abs(B_opt)>10^(-3);
        
        [FP3,TP3] = evaluate_error_connectivity(B_opt);
        FP_tot = FP_tot+FP; TP_tot = TP_tot+TP;
        FP_tot_deconfound = FP_tot_deconfound+FP3; TP_tot_deconfound = TP_tot_deconfound+TP3;
        
    end
    FP_final(:,j) = FP_tot/10; TP_final(:,j) = TP_tot/10;
    
    FP_final_deconfound(:,j) = FP_tot_deconfound/10; TP_final_deconfound(:,j) = TP_tot_deconfound/10;
    
end