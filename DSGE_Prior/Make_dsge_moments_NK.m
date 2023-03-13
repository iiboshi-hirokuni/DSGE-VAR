function [XXdum, XYdum, YYdum] = Make_dsge_moments_NK(para,ZZ,p)

% This procedure evaluates the likelihood function of the 
% monetary DSGE model
% retcode = -1 : non existence
%         = 0  : existence and uniqueness
%         = 1  : existence and non-uniqueness

% Parameters
%  h,  sigma_c, sigma_L, beta, phi, tau, Rk, gam_p, 
%  gam_w, psi_p, psi_w, alpha, psi, k_y, g_y,
%  rho_m, mu_pi, mu_y
%  e_c, e_inv, e_q, e_L, e_w, e_z, e_p, e_g, e_m

% npara   = size(para, 1);
retcode = 0;

% solve the DSGE model

[T1,TC,TEPS,RC] = dsgesolv_NK(para);
%DD = zeros(size(ZZ, 1),1);

nseries  = size(ZZ, 1);
nstate   = size(T1, 2);
nshock  = size(TEPS,2);

% nobs      = size(yy, 1);
% loglh     = 0;
% loglhzero = -1E8;
% obsmean   = zeros(nobs, nseries);
% obsvar    = zeros(nobs, nseries);
% nu_save   = zeros(nseries,nobs);
% ft_save   = zeros(nseries,nseries*nobs); 
% kg_save   = zeros(nstate,nseries*nobs);
% pt_save   = zeros(nstate,nstate*nobs);
% shock     = zeros(nobs, nshock);

% /**********************************************************
% ** Check determinacy 
% **********************************************************/
if (RC(1) == 1) && (RC(2)==1);
   %/* determinacy */
   retcode(1) = 0;
   TT = T1;
   RR = TEPS;
   
elseif (RC(1) == 1) && (RC(2)==0) 
   %/* indeterminacy */
   retcode(1) = 1;
   TT = T1;
   RR = TEPS;
   rloglh = loglhzero;
   return;

else
   %/* no equilibrium exists, numerical problems */
   retcode(1) = RC(1);
   rloglh = loglhzero;
   return;

end

% create system matrices for state space model

% These matrices are regime independent

DD = zeros(nseries,1);

z_star_bar = para(14,1);
psi_bar    = para(15,1);

% l_bar   = mean(yy(:,5));  %1.0;
% r_n_bar = mean(yy(:,8));
% r_l_bar = mean(yy(:,10));
pi_bar  = 0.25;
r_n_bar  = 1;

DD = zeros(nseries,1);
DD(1, 1) = z_star_bar;
DD(2, 1) = z_star_bar;
DD(3, 1) = z_star_bar+psi_bar;
DD(4, 1) = z_star_bar;
DD(5, 1) = pi_bar;
DD(6, 1) = r_n_bar;


% HH = zeros(nseries,nseries);
HH = 0.01*eye(nseries);
QQ = createcov(para(31:39,1))/100;
% VV = zeros(nshock,nseries);

% Check whether covariance matrix QQ is positive definite

% if sum(eig(QQ) < 0) > 0
%    loglh = loglhzero;
%    return;
% end

% We can now define the initial mean and variance for the state vector
%
% At = zeros(nstate*2,1);
% At(1)=yy(1,1); At(2)=yy(1,2); At(3)=yy(1,3); At(4)=yy(1,4);    
% At(5)=yy(1,5); At(6)=yy(1,6); At(7)= yy(1,7); 
% At(8)=yy(1,8); At(9)=yy(1,9); At(10)= yy(1,10); 


Pt = dlyap(TT,RR*QQ*RR');
Pt = [Pt zeros(nstate);
      zeros(nstate) Pt];

TT = [TT zeros(nstate);
      diag(ones(nstate,1)) zeros(nstate)];

RR = [RR; zeros(nstate,nshock)];

if p == 1
  YYdum = ZZ*Pt*ZZ'+ HH + DD*DD';
  XYdum = ZZ*TT*Pt*ZZ'+ DD*DD';
  XYdum = [XYdum; DD'];
  XXdum = [ YYdum DD; DD' 1];
elseif p == 2  
  YYdum = ZZ*Pt*ZZ'+ HH + DD*DD';
  XX1dum = ZZ*TT*Pt*ZZ'+ DD*DD';
  XX2dum = ZZ*TT*TT*Pt*ZZ'+ DD*DD'; 
  
  XYdum = [XX1dum; XX2dum; DD'];
  
  XXdum = [ YYdum XX1dum DD;...
            XX1dum YYdum DD; ...            
            DD'     DD'   1];   
  
elseif p == 4  
  YYdum = ZZ*Pt*ZZ'+ HH + DD*DD';
  XX1dum = ZZ*TT*Pt*ZZ'+ DD*DD';
  XX2dum = ZZ*TT*TT*Pt*ZZ'+ DD*DD';
  XX3dum = ZZ*TT*TT*TT*Pt*ZZ'+ DD*DD';
  XX4dum = ZZ*TT*TT*TT*TT*Pt*ZZ'+ DD*DD';
  
  XYdum = [XX1dum; XX2dum; XX3dum; XX4dum; DD'];
  
  XXdum = [ YYdum XX1dum XX2dum XX3dum DD;...
            XX1dum YYdum XX1dum XX2dum DD; ...
            XX2dum XX1dum YYdum XX1dum DD; ...
            XX3dum XX2dum XX1dum YYdum DD;...
            DD'    DD'    DD'    DD'   1];
end


% 
% %-----------------------------
function [omega] = createcov(para)

npara = max(size(para));
omega = zeros(npara, npara);
for i = 1:npara
  omega(i, i) = para(i)^2;
end
