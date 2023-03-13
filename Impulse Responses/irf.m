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
 
% %  �T�C������  �̒�`
% % if y(2,1:j)<=r & y(3,1:j)>=r & y(4,1:j)<=r;
% %if  y(3,1:j)>=r & y(4,1:j)<=r;    
%     % ���Z����V���b�N�ɑ΂��ăC���t���� y(2,) �̔����� j���܂Ń}�C�i�X
%     % ���Z����V���b�N�ɑ΂��ċ��� y(3,) �̔����� j���܂Ńv���X
%     % ���Z����V���b�N�ɑ΂��Ď����ݕ� y(4,) �̔����� j���܂Ń}�C�i�X
%     
%     d=1;
% else
%     d=0;
% end
% end
% end
