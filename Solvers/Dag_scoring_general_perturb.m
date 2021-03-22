function[estim_params,likelihood_score] = Dag_scoring_general_perturb(Dag_adj)

global Sigma
global nvec
global learn_rate
global perturbed_latent
global latent
global psi_max
global rank_est
global R
global ratio_BIC


% first initialize connectivity matrix
p = size(Dag_adj,1);
B_init = zeros(p);
for i = 1:p
    ind = find(Dag_adj(i,:) ~= 0);
    B_init(i,ind) = Sigma{1}(i,ind)*(Sigma{1}(ind,ind))^(-1);
end

% initialize Gamma_init, do_init
[U,D,~] = svd((eye(p)-B_init)*Sigma{1}*(eye(p)-B_init)');
%[U D V] = svd(randn(p));
Gamma_init = latent*U(:,1:rank_est)*D(1:rank_est,1:rank_est)^(1/2);
%d0_init = diag(U(:,rank_est+1:end)*D(1+rank_est:end,1+rank_est:end)*...
 %   U(:,rank_est+1:end)');
d0_init = diag((eye(p)-B_init)*Sigma{1}*(eye(p)-B_init)');

% initialize psi,zeta
psi_init = zeros(length(nvec),1);
zeta_init = zeros(length(nvec),1);
for j = 1:length(nvec)
    psi_vec = linspace(0.01,psi_max,200);
    zeta_vec = linspace(0.01,10,200);
    h = 1;
    for k = 1:length(psi_vec)
        for l = 1:length(zeta_vec)
            M = (Gamma_init*Gamma_init')*(1+psi_vec(k))+diag(d0_init)+...
                zeta_vec(l)*eye(p);
            score(h) = log(det(M))+trace(M^(-1)*((eye(p)-diag(R(:,j))*B_init)*Sigma{j}*(eye(p)-diag(R(:,j))*B_init)'));
            param_vec(h,1:2) = [k,l];
            h = h+1;
        end
    end
    [~,ind] = min(score);
    psi_init(j) = psi_vec(param_vec(ind,1));
    zeta_init(j) = zeta_vec(param_vec(ind,2));
end
psi_init = psi_init*perturbed_latent;

%% Perform alternating B and the rest (solved via gradient descent)
% initialize parameters
%psi_init(2:end) = 1+rand(length(nvec)-1,1);
B_new = B_init; B_old = rand(p,p);
d0_new = d0_init; d0_old = rand(p,1);
Gamma_new = Gamma_init; Gamma_old = rand(p,rank_est);
psi_init(1) = 0;
psi_new = psi_init; psi_old = rand(length(nvec),1);

for j = 2:length(nvec)
    zeta_new(1:p,j) = zeta_init(j)*ones(p,1);
end

if length(nvec) < 2
    zeta_new = rand(p,length(nvec));
end
zeta_old = rand(p,length(nvec));


% max_dev checks closeness of successive parameters
max_dev = 100;
fluctuation_count = 0;
sparse_pattern = find(abs(Dag_adj)<= 10^(-3));
while max_dev >= 10^(-2) 
    
    max_dev_old = max_dev;
    B_old = B_new; Gamma_old = Gamma_new; d0_old = d0_new;
    psi_old = psi_new; zeta_old = zeta_new;
    
    M = cell(length(nvec),1);
    for j = 1:length(nvec)
        M{j} = (diag(R(:,j))*((Gamma_old*Gamma_old')*(1+psi_old(j))+diag(d0_old)+...
            diag(zeta_old(:,j)))*diag(R(:,j)) + diag(1-R(:,j))*((Gamma_old*Gamma_old')*(1+psi_old(j))+diag(d0_old)+...
            diag(zeta_old(:,j)))*diag(1-R(:,j)))^(-1);
    end
    
    B_new = full(solve_B_fixed_everything(sparse_pattern,M));
    
    
    %% perform gradient descent for a fixed B
    overall_grad = 0.1;
    Gamma_current = Gamma_old;
    d0_current = d0_old;
    zeta_current = zeta_old;
    psi_current = psi_old;
    score = [];
    increase_lik_sc = 0;
     iter = 0;
     score_improve = 1;

    while (score_improve >= 10^(-6) )
        
        iter = iter+1;
        Gradient_Gamma = zeros(p,rank_est);
        Gradient_L = zeros(p,p);
        Gradient_d0 = zeros(p,1);
        Gradient_zeta = zeros(p,length(nvec));
        Gradient_psi = zeros(length(nvec),1);
        
        for j = 1:length(nvec)
            
            % set parameters that are used often
            M =  (diag(R(:,j))*((1+psi_current(j))*(Gamma_current*Gamma_current')...
                +diag(d0_current)+(diag(zeta_current(:,j))))*diag(R(:,j))+...
                diag(1-R(:,j))*((1+psi_current(j))*(Gamma_current*Gamma_current')...
                +diag(d0_current)+(diag(zeta_current(:,j))))*diag(1-R(:,j)))^(-1);
            N = (eye(p)-diag(R(:,j))*B_new)*Sigma{j}*(eye(p)-diag(R(:,j))*B_new)';

            
            % Ensure positive definite matrix
            if min(eig(M))<0
                break;
            end
            
            
            % computing gradient for Gamma
            for k = 1:p
               if R(k,j) == 1
                for l = 1:rank_est
                    temp = zeros(p,p);
                    temp(k,k) = Gamma_current(k,l)'*(1+psi_current(j));
                    temp(k,setdiff(1:p,k)) = Gamma_current(setdiff(1:p,k),l)'*(1+psi_current(j));
                    temp = temp+temp';
                    Gradient_Gamma(k,l) = Gradient_Gamma(k,l)+...
                        nvec(j)/nvec(1)*trace(M*temp)-nvec(j)/nvec(1)*trace(M*N*M*temp);
                end
              end
            end
            
            % computing gradient for d0
            Gradient_d0 = Gradient_d0+nvec(j)/nvec(1)*diag(M-(M*N*M));
            
            % computing gradient for zeta,psi
            H = diag(M-(M*N*M));
            if j > 1
                Gradient_zeta(R(j),j) = Gradient_zeta(R(j),j)+nvec(j)/nvec(1)*diag(H(R(j),R(j)));
                Gradient_psi(j) = Gradient_psi(j) + nvec(j)/nvec(1)*(trace(M*(Gamma_current*Gamma_current'))...
                    -trace(N*M*(Gamma_current*Gamma_current')*M));
            end
            ind_nondo = find(R(:,j) == 1);
            score_diag(j) = nvec(j)/nvec(1)*(-log(det(M(ind_nondo,ind_nondo)))+...
                trace(M(ind_nondo,ind_nondo)*N(ind_nondo,ind_nondo)));

        end
        
        
        %update parameters
        Gamma_current = latent*(Gamma_current - learn_rate* Gradient_Gamma);
        d0_current = (d0_current - learn_rate*Gradient_d0).*((d0_current>=0));
        zeta_current = zeta_current - learn_rate*Gradient_zeta; zeta_current = zeta_current.*(zeta_current>=0);
        psi_current = psi_current - learn_rate*Gradient_psi;psi_current = psi_current.*(psi_current>=0);
        psi_current(find(psi_current>=psi_max))=psi_max;
        psi_current = psi_current*perturbed_latent;
        
        %diagnostics
          score = [score sum(score_diag)];
            if iter >1 
                score_diff = (score(end)-score(end-1));
                if score_diff > 10^(-5)
                    increase_lik_sc = increase_lik_sc+1;
                    %error('Cannot calculate with given values')
                end
                score_improve = abs((score(end)-score(end-1))/score(end-1));
            end
            [norm(Gradient_Gamma), norm(Gradient_d0),norm(Gradient_zeta),...
            norm(Gradient_psi)];
        
        %compute norm of gradients for stopping criteria
        overall_grad = max([norm(Gradient_Gamma), norm(Gradient_d0),norm(Gradient_zeta),norm(Gradient_psi)]);
        
    end
    
    % update parameters
    Gamma_new = Gamma_current; d0_new = d0_current; zeta_new = zeta_current; psi_new = psi_current;
    
    % compute closeness of successive estimates
    max_dev = max([max(abs(B_old-B_new))]);
    if max_dev_old < max_dev
        fluctuation_count = fluctuation_count+1;
    end
end

estim_params = cell(5,1);estim_params{1}=B_new;estim_params{2}=Gamma_new;
estim_params{3} = d0_new; estim_params{4} = zeta_new; estim_params{5} = psi_new;

S = (eye(p)-B_new)'*(eye(p)-B_new);
S = S-diag(diag(S));


% compute likelihood score
likelihood_score = 0;
for j = 1:length(nvec)
    M = (diag(R(:,j))*((1+psi_new(j))*(Gamma_new*Gamma_new')...
                +diag(d0_new)+diag(zeta_new(:,j)))*diag(R(:,j))+...
                diag(1-R(:,j))*((1+psi_new(j))*(Gamma_new*Gamma_new')...
                +diag(d0_new)+diag(zeta_new(:,j)))*diag(1-R(:,j)))^(-1);
    N = (eye(p)-diag(R(:,j))*B_new)*Sigma{j}*(eye(p)-diag(R(:,j))*B_new)';
    ind_nondo = find(R(:,j) == 1);
    likelihood_score = likelihood_score+ nvec(j)/sum(nvec)*(-log(det(M(ind_nondo,ind_nondo)))+...
        trace(M(ind_nondo,ind_nondo)*N(ind_nondo,ind_nondo)))+ratio_BIC*log(sum(nvec))/(sum(nvec))*...
        length(find(abs(S)>=10^(-3)))/2;
end