

close all
clear all
clc

p          = 2;                 % Number of lags in the VAR

lambda_save = [0.001, 0.005 0.1 0.2  0.3 0.35 0.4 0.45 0.5 0.75 1 2  ]';
mu_save = [0 0.2 0.4 0.5 0.6 0.8  1]';

mdd_save = zeros(size(lambda_save,1),size(mu_save,1));

for itr_lambda = 1:size(lambda_save,1)
    for itr_mu = 1:size(mu_save,1)

        lambda = lambda_save(itr_lambda)
        ratio =  lambda/(1+lambda)
        
        mu     =  mu_save(itr_mu)
        
        DSGE_VAR

        mdd_save(itr_lambda,itr_mu)=mdd;
    end
end

ratio = zeros(size(lambda_save,1),1);

for i = 1:size(lambda_save,1)
 ratio(i) = lambda_save(i)/(1+lambda_save(i));
end
 
save 'mdd_save.mat' 'mdd_save' 'mu_save' 'ratio';

Plot_Graph;

