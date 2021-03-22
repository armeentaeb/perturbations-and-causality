function[B_opt] = learning_and_prediction(equivalent_DAGs,Test_Data)
global ratio_BIC
global nvec

ratio_BIC = 0;
[optimal_Bv1,Gamma,d0,scorev1] = score_all_dags_set(equivalent_DAGs);

[sortval ind1] = sort(scorev1,'ascend');
[ind2] = find(abs((sortval-sortval(1))/sortval(1))<=10^(-3));

for opt_iter = 1:length(ind2)


%[val1,ind1] = min(scorev1);
M = optimal_Bv1{ind1(ind2(opt_iter))};
%M = optimal_Bv1{ind1};
[row col] = find(abs(optimal_Bv1{ind1(ind2(opt_iter))})>0.001);
vec = zeros(length(row),1);
for i = 1:length(row)
    vec(i) = M(row(i),col(i));
end
[val ind] = sort(abs(vec),'ascend');
row = row(ind);
col = col(ind);


    for l = 1:length(row)
        clear equivalent_DAGs2;
        equivalent_DAGs2 = [];
        M(row(l),col(l))= 0;
        equivalent_DAGs2{1} = M;
        ratio_BIC = 0;
        [optimal_Bv2,Gamma,d0,scorev2] = score_all_dags_set(equivalent_DAGs2');
        p = size(optimal_Bv2{1},1);
        train(l,opt_iter) = scorev2;
        test_sigma = cov(Test_Data);
        precision= (eye(p)-optimal_Bv2{1})'*(diag(d0{1})+Gamma{1}*Gamma{1}')^(-1)*(eye(p)-optimal_Bv2{1});
        test(l,opt_iter) = -log(det(precision))+trace(precision*test_sigma);
          M = optimal_Bv2{1};
        score_val(l,opt_iter) = scorev2;
        DAG{l,opt_iter} =  M;
        result = M;
        score1 = scorev2;
    end

end    
    
    
    
DAG_final = cell(length(BIC_thresh),1);
for k = 1:length(BIC_thresh)  
    [val ind] = min(train(k,:));
    test_sel(k) = test(k,ind);
    DAG_final{k} = DAG{ind};
end

[val ind] = min(test_sel);
B_opt = DAG_final{ind};


    