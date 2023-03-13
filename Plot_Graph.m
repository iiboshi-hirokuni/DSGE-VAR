
clear all
clc 

load 'mdd_save.mat';

figure(4000)
 [X, Y] = meshgrid(mu_save, ratio);
 
%  mesh(X,Y,mdd_save);
z =mdd_save;
      
 surface(X,Y,mdd_save,gradient(z)); 
  rotate3d;
   ylabel('\lambda/(1+\lambda)');
   xlabel('\mu');
   zlabel('Marginal Likelihood')
  view(30+20, 20);
  title('Posteror Distribution of Hyperparameters')


 figure(5000)
 hold on
 for i = 1:size(mu_save,1)      
     plot(ratio(3:11), mdd_save(3:11,i)','LineWidth',1.5)
 end 
 hold off
   xlabel('\lambda/(1+\lambda)')
   ylabel('Marginal Likelihood')
    ylim([-380 -340])
    xlim([0.15 0.4])
   legend('\mu=0 (FF)', '\mu=0.2', '\mu=0.4', '\mu=0.5', '\mu=0.6', '\mu=0.8', '\mu=1(NK)')