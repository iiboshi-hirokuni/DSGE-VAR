function [y] = irf(Phi,Sigma,H,j)

% This function draw the IRFs from the posterior density. The matrix Omega
% is drawn from the unit ball untill the sign restrictions on the impulse
% response of Inflation(y(2,:)), the Interest Rate (y(3,:)) and Real Money 
% (y(4,:)) are satisfied for j periods. 

[m,n] = size(Sigma);
y = ones(n,H);
d=0;
r = zeros(1,j);

% while d==0
%e = randn(n,1);
%  i.e., Omega * Omega' = I
%Omega = (e./norm(e));

%  Omega = QR decompotion
% RWZ (2010) p688, algorithm 2 
% x = randn(n,n);
% [P, R]=qr(x);    %  i.e., Q * Q' = I
% Omega = P*shock;
shock =[1;0;0;0;0;0];
Omega = shock;


for i=1:H
    y(:,i) = impulseresponce(Phi,Sigma,Omega,i);
end
 
% %  サイン制約  の定義
% % if y(2,1:j)<=r & y(3,1:j)>=r & y(4,1:j)<=r;
% %if  y(3,1:j)>=r & y(4,1:j)<=r;    
%     % 金融政策ショックに対してインフレ率 y(2,) の反応は j期までマイナス
%     % 金融政策ショックに対して金利 y(3,) の反応は j期までプラス
%     % 金融政策ショックに対して実質貨幣 y(4,) の反応は j期までマイナス
%     
%     d=1;
% else
%     d=0;
% end
% end
% end
