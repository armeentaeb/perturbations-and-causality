function [B] = solve_B_fixed_everything(sparse_pattern,M)
global Sigma
global nvec
global R;
p = size(M{1},1);
cvx_begin quiet
variable B(p,p)
variable t
minimize t
subject to 
v = 0;
for j = 1:length(nvec)
    ind_do = find(R(:,j) == 1); F = eye(p); F = F(ind_do,:);
    v = sum(sum_square_abs(M{j}^(1/2)*(eye(p)-diag(R(:,j))*B)*Sigma{j}^(1/2)))+v;
end
v <= t
B(sparse_pattern) == 0;
cvx_end


%% solve using projected gradient descent

% Bnew = rand(p,p);Bnew(sparse_pattern) = 0;
% grad = 1;
% learn_rate = 0.1;
% while norm(grad)>10^(-5)
%    
%     grad = zeros(p);
%     for j = 1:length(nvec)
%     grad = grad + 2*nvec(j)/nvec(1)*M{j}*(Bnew-eye(p))*Sigma{j};
%     end
%     Bnew = Bnew - learn_rate*grad;
%     Bnew(sparse_pattern) = 0;
%     norm(grad)
% end








