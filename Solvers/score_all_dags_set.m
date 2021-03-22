function[Bfinal,Gamma,d0,likelihood_score] = score_all_dags_set(DAG_set)
p = size(DAG_set{1},1);


score_init_v1 = 1000000;
score_init_v2 = 1000000;

for j = 1:size(DAG_set,1)
    tic
    [estimated_param, likelihood] = Dag_scoring_general_perturb(DAG_set{j});
    toc
    likelihood_score(j) = likelihood;
    Bfinal{j} = estimated_param{1};
    Gamma{j} = estimated_param{2};
    d0{j} = estimated_param{3};
end

