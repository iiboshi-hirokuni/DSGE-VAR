function [y] = irf(Phi,Sigma,H,j,P_FF,P_NK,mu,shock_no)

% This function draw the IRFs from the posterior density. The matrix Omega
% is drawn from the unit ball untill the sign restrictions on the impulse
% response of Inflation(y(2,:)), the Interest Rate (y(3,:)) and Real Money 
% (y(4,:)) are satisfied for j periods. 

[m,n] = size(Sigma);
y = zeros(n,H);
d=0;
r = zeros(1,j);

% while d==0
% %e = randn(n,1);
% %  i.e., Omega * Omega' = I
% %Omega = (e./norm(e));

%  Omega = QR decompotion
% RWZ (2010) p688, algorithm 2 
% x = randn(n,n);
% [P, R]=qr(x);    %  i.e., Q * Q' = I
if shock_no ==1  %preference
   shock =[-1;0;0;0;0;0];
elseif shock_no ==6  % monetary policy 
  shock =[0;0;0;0;0;1];
end

Omega = (mu*P_NK+(1-mu)*P_FF)*shock;
% Sigma=eye(6);

for i=1:H
    dy(:,i) = impulseresponce(Phi,Sigma,Omega,i);
 
    if i>1
%       y(1:4,i) = y(1:4,i-1)+dy(1:4,i);
      y(1:4,i) = dy(1:4,i);
     y(5:end,i) = dy(5:end,i);
    else
     y(:,1)=dy(:,1);   
    end
end
 

