%
%  Calculationg of inefficiency
% 
%
%

sum_rho =0;

num = nsim-nburn;
B_m = round(num/100);   % Bandwidth of Windows  


inefficiency= zeros(n*p+1,1);

for m = 1:1:n*p+1
 for j=1:1:B_m

  %  Parzen window (See Kim, Shaphard and Chib (1998) )  
    z = j/B_m;
   if (z < 0.5) 
      K = 1 - 6*z^2+6*z^3; 
   elseif (z >= 0.5)
      K=2*(1-z)^3;
   end    
      
    rho_hat = corrcoef( [para_save(1+j:num,m) para_save(1:num-j,m)] );

    sum_rho = sum_rho + K*rho_hat(2,1);
    
 end

inefficiency(m) = 1 + 2*B_m/(B_m-1)*sum_rho; 

end

