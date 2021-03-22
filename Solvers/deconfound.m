function[precision,graphModelLVGM] = deconfound(X)
temp = pwd;
root = temp(1:36);
%root = '/Users/armeentaeb/Dropbox/Research_at_ETH/Code';
addpath(strcat(root,'/Solvers/logdetppa-0'))
addpath(strcat(root,'/Solvers/logdetppa-0/solver'))
addpath(strcat(root,'/Solvers/logdetppa-0/util'))
addpath(strcat(root,'/Solvers/logdetppa-0/solver/mexfun'))
global LVGM; % Latent variable graphical modeling
LVGM = 2;
global R;
global Sigma
global rank_est
% 

%train and test amount
n_train = size(X,1)*8/10;
n_test = size(X,1)-n_train;
p = size(X,2);
X_train = X;

global TrainTest;
TrainTest = cell(4,1);
TrainTest{1} = X(1:n_train,:);
TrainTest{2} = X(n_train+1:end,:);
TrainTest{3} = [];
TrainTest{4} = []; 
num_bags = 20;
alpha = 0.7;

lambdaVec = linspace(0.1,0.05,30);
gammaVec = linspace(0.5,0.01,30);

h = 1;
for i = 1:length(lambdaVec)
    for j = 1:length(gammaVec)
        [precision, lowRankLVGM1, graphModelLVGM1] = ObtainEstimate(TrainTest,LVGM,[lambdaVec(i), gammaVec(j)]);
        params(h,1:2) = [lambdaVec(i), gammaVec(j)];
        sparse(h) = sum(sum(abs(graphModelLVGM1)>=10^(-3)));
        test_error(h) = -log(det(precision))+trace(precision*cov(TrainTest{2}));
        rank_comp(h) = length(find(svd(lowRankLVGM1)>=10^(-3)));
        h = h+1;
    end
end
ind = find(rank_comp <= rank_est);
[~,ind2] = min(test_error(ind));
[precision, lowRankLVGM1, graphModelLVGM] = ObtainEstimate(TrainTest,LVGM,[params(ind(ind2(1)),1), params(ind(ind2(1)),2)]);
















