function [T1,TC,T0, RC] = dsgesolv_SW(para)
%^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   Solve New Keynesian Business Cycle Model
%^^^^^^^^^^^^^^^^^^^^^^^^^^^

% assign names to parameters
sigma      = para(1,1);
theta      = para(2,1);
kai        = para(3,1);
inv_zeta   = para(4,1);
mu         = para(5,1);
phi_o_y    = para(6,1);
gamma_w    = para(7,1);
ksi_w      = para(8,1);
gamma_p    = para(9,1);
ksi_p      = para(10,1);
phi_r      = para(11,1);
phi_pi     = para(12,1);
phi_y      = para(13,1);
z_star_bar = para(14,1);
psi_bar    = para(15,1);
%eta        = para(16,1);
% n_k        = para(17,1);
% mu_E       = para(18,1);
% % r_E_bar    = para(19,1);

rho_b      = para(20,1);
rho_g      = para(21,1);
rho_w      = para(22,1);
rho_p      = para(23,1);
rho_r      = para(24,1);
rho_nu     = para(25,1);
rho_z      = para(26,1);
rho_psi    = para(27,1);
rho_i      = para(28,1);
% rho_efp    = para(29,1);
% rho_nw     = para(30,1);

% sigma_b    = para(31,1);
% sigma_g    = para(32,1);
% sigma_w    = para(33,1);
% sigma_p    = para(34,1);
% sigma_r    = para(35,1);
%  sigma_nu   = para(36,1);
% sigma_z    = para(37,1);
% sigma_psi  = para(38,1);
% sigma_i  =   para(39,1);
% sigma_efp  = para(40,1);
% sigma_nw   = para(41,1);

delta    = 0.06;
alpha    = 0.37;
lambda_w = 0.2;
lambda_i = 0.2;

c_y    = 1 - alpha;
i_y    = alpha;
% l_bar = 1.0;
r_n_bar= 1.0;
pi_bar = 0.25;

z_star = exp(0.01*z_star_bar);
psi    = exp(0.01*psi_bar);
% l = exp(0.01*l_bar);
pii    = exp(0.01*pi_bar);
r_n    = exp(0.01*r_n_bar);
% r_E    = exp(0.01*r_E_bar);

% define matrices of canonical system

neq  = 35;      %  Num of stable and unstable Variables
nex  = 9;       %  Num of Shock
nend = 9;       %  Num of Unstable Variables

GAM0j = zeros(neq,neq);
GAM1j = zeros(neq,neq);
C = zeros(neq,1);
PSI0j = zeros(neq,nex);
PPIj = zeros(neq,nend);

  GAM0j(1,1) = sigma/((theta/z_star-1.0)*((pii*theta)/r_n- ...
     1.0))+(pii*sigma*theta^2)/(r_n*z_star*(theta/z_star-1.0 ...
     )*((pii*theta)/r_n-1.0));
      GAM0j(1,2) = 1.0;
      GAM0j(1,15) = 1.0/((pii*theta)/r_n-1.0);
      GAM0j(1,23) = (sigma*theta)/(z_star*(theta/z_star-1.0)*(... 
     (pii*theta)/r_n-1.0));
      GAM0j(1,28) = -(pii*sigma*theta)/(r_n*(theta/z_star-1.0) ...
     *((pii*theta)/r_n-1.0));
      GAM0j(1,34) = -(pii*theta)/(r_n*((pii*theta)/r_n-1.0));
      GAM0j(1,35) = -(pii*sigma*theta)/(r_n*(theta/z_star-1.0) ...
      *((pii*theta)/r_n-1.0));
      GAM0j(2,2) = 1.0;
      GAM0j(2,14) = -1.0;
      GAM0j(2,29) = -1.0;
      GAM0j(2,32) = 1.0;
      GAM0j(2,35) = sigma;
      GAM0j(3,2) = ((ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))/(ksi_w ...
      *((kai*(lambda_w+1.0))/lambda_w+1.0));
      GAM0j(3,3) = (pii*z_star)/r_n+((ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))...
          /(ksi_w*((kai*(lambda_w+1.0))/lambda_w+1.0))+1.0;
      GAM0j(3,8) = (gamma_w*pii*z_star)/r_n+1.0;
      GAM0j(3,15) = -((ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))/(ksi_w ...
          *((kai*(lambda_w+1.0))/lambda_w+1.0));
      GAM0j(3,17) = -1.0;
      GAM0j(3,23) = 1.0;
      GAM0j(3,24) = -(kai*(ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))/ ...
         (ksi_w*((kai*(lambda_w+1.0))/lambda_w+1.0));
      GAM0j(3,30) = -(pii*z_star)/r_n;
      GAM0j(3,32) = -(pii*z_star)/r_n;
      GAM0j(3,35) = -(pii*z_star)/r_n;
      GAM0j(4,4) = (delta-1.0)/(psi*r_n);
      GAM0j(4,5) = -(delta-1.0)/(psi*r_n)-1.0;
      GAM0j(4,22) = 1.0;
      GAM0j(4,26) = 1.0;
      GAM0j(5,3) = -1.0;
      GAM0j(5,5) = 1.0;
      GAM0j(5,6) = 1.0;
      GAM0j(5,22) = -1.0;
      GAM0j(5,23) = -1.0;
      GAM0j(5,24) = -1.0;
      GAM0j(6,4) = mu;
      GAM0j(6,5) = -mu;
      GAM0j(6,6) = 1.0;
      GAM0j(7,3) = alpha-1.0;
      GAM0j(7,5) = -alpha;
      GAM0j(7,7) = 1.0;
      GAM0j(8,7) = -((ksi_p-1.0)*((ksi_p*pii*z_star)/r_n-1.0))/ksi_p;
      GAM0j(8,8) = (gamma_p*pii*z_star)/r_n+1.0;
      GAM0j(8,18) = -1.0;
      GAM0j(8,32) = -(pii*z_star)/r_n;
      GAM0j(9,6) = -alpha*(phi_o_y+1.0);
      GAM0j(9,11) = 1.0;
      GAM0j(9,22) = alpha*(phi_o_y+1.0);
      GAM0j(9,23) = alpha*(phi_o_y+1.0);
      GAM0j(9,24) = (alpha-1.0)*(phi_o_y+1.0);
      GAM0j(10,1) = -c_y;
      GAM0j(10,11) = 1.0;
      GAM0j(10,12) = -i_y;
      GAM0j(10,16) = -1.0;
      GAM0j(11,6) = (delta+(psi*r_n)/pii-1.0)/(psi*z_star);
      GAM0j(11,12) = -(delta-1.0)/(psi*z_star)-1.0;
      GAM0j(11,13) = 1.0;
      GAM0j(11,19) = -(delta-1.0)/(psi*z_star)-1.0;
      GAM0j(11,22) = -(delta-1.0)/(psi*z_star);
      GAM0j(11,23) = -(delta-1.0)/(psi*z_star);
      GAM0j(12,4) = 1.0;
      GAM0j(12,12) = -inv_zeta-(inv_zeta*pii*z_star)/r_n;
      GAM0j(12,19) = 1.0;
      GAM0j(12,22) = -inv_zeta;
      GAM0j(12,23) = -inv_zeta;
      GAM0j(12,27) = -1.0;
      GAM0j(12,33) = (inv_zeta*pii*z_star)/r_n;
      GAM0j(12,35) = (inv_zeta*pii*z_star)/r_n;
      GAM0j(13,8) = phi_pi*(phi_r-1.0)*(1.0/4.0);
      GAM0j(13,9) = phi_pi*(phi_r-1.0)*(1.0/4.0);
      GAM0j(13,10) = phi_pi*(phi_r-1.0)*(1.0/4.0);
      GAM0j(13,11) = phi_y*(phi_r-1.0);
      GAM0j(13,14) = 1.0;
      GAM0j(13,20) = -1.0;
      GAM0j(14,14) = -1.0;
      GAM0j(14,25) = 1.0;
      GAM0j(14,32) = 1.0;
      GAM0j(15,1) = 1.0;
      GAM0j(16,2) = 1.0;
      GAM0j(17,3) = 1.0;
      GAM0j(18,26) = 1.0;
      GAM0j(19,5) = 1.0;
      GAM0j(20,8) = 1.0;
      GAM0j(21,12) = 1.0;
      GAM0j(22,15) = 1.0;
      GAM0j(23,23) = 1.0;
      GAM0j(24,15) = 1.0;
      GAM0j(25,16) = 1.0;
      GAM0j(26,17) = 1.0;
      GAM0j(27,18) = 1.0;
      GAM0j(28,19) = 1.0;
      GAM0j(29,20) = 1.0;
      GAM0j(30,21) = 1.0;
      GAM0j(31,22) = 1.0;
      GAM0j(32,27) = 1.0;
      GAM0j(33,21) = -1.0;
      GAM0j(33,22) = alpha/(alpha-1.0);
      GAM0j(33,23) = 1.0;
      GAM0j(34,9) = 1.0;
      GAM0j(35,10) = 1.0;

   %%
    GAM1j(1,1) = (sigma*theta)/(z_star*(theta/z_star-1.0)*(( ...
      pii*theta)/r_n-1.0));
      GAM1j(3,3) = 1.0;
      GAM1j(3,8) = gamma_w;
      GAM1j(4,4) = -1.0;
      GAM1j(5,13) = -1.0;
      GAM1j(8,8) = gamma_p;
      GAM1j(9,13) = alpha*(phi_o_y+1.0);
      GAM1j(11,13) = -(delta-1.0)/(psi*z_star);
      GAM1j(12,12) = -inv_zeta;
      GAM1j(13,10) = phi_pi*(phi_r-1.0)*(-1.0/4.0);
      GAM1j(13,11) = phi_y*(phi_r-1.0);
      GAM1j(13,14) = phi_r;
      GAM1j(15,28) = 1.0;
      GAM1j(16,29) = 1.0;
      GAM1j(17,30) = 1.0;
      GAM1j(18,25) = 1.0;
      GAM1j(19,31) = 1.0;
      GAM1j(20,32) = 1.0;
      GAM1j(21,33) = 1.0;
      GAM1j(22,34) = 1.0;
      GAM1j(23,35) = 1.0;
      GAM1j(24,15) = rho_b;
      GAM1j(25,16) = rho_g;
      GAM1j(26,17) = rho_w;
      GAM1j(27,18) = rho_p;
      GAM1j(28,19) = rho_nu;
      GAM1j(29,20) = rho_r;
      GAM1j(30,21) = rho_z;
      GAM1j(31,22) = rho_psi;
      GAM1j(32,27) = rho_i;
      GAM1j(34,8) = 1.0;
      GAM1j(35,9) = 1.0;

      %%
     PSI0j(24,1) = 1.0;
      PSI0j(25,2) = 1.0;
      PSI0j(26,3) = 1.0;
      PSI0j(27,4) = 1.0;
      PSI0j(28,5) = 1.0;
      PSI0j(29,6) = 1.0;
      PSI0j(30,7) = 1.0;
      PSI0j(31,8) = 1.0;
      PSI0j(32,9) = 1.0;

      %%
      PPIj(15,1) = 1.0;
      PPIj(16,2) = 1.0;
      PPIj(17,3) = 1.0;
      PPIj(18,9) = 1.0;
      PPIj(19,4) = 1.0;
      PPIj(20,5) = 1.0;
      PPIj(21,6) = 1.0;
      PPIj(22,7) = 1.0;
      PPIj(23,8) = 1.0;

% QZ(generalized Schur) decomposition by GENSYS
%[T1,TC,T0,TY,M,TZ,TETA,GEV,RC) = gensys(GAM0,GAM1,C,PSI0,PPI,1,1);
[T1,TC,T0,fmat,fwt,ywt,gev,RC,loose] = gensys(GAM0j,GAM1j,C,PSI0j,PPIj,1);
