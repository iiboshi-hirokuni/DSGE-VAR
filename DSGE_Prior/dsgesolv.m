function [T1,TC,T0, RC] = dsgesolv(para)
%******************************************************
%   Solve New Keynesian Business Cycle Model
%******************************************************

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
eta        = para(16,1);
 n_k        = para(17,1);
 mu_E       = para(18,1);
 r_E_bar    = para(19,1);

rho_b      = para(20,1);
rho_g      = para(21,1);
rho_w      = para(22,1);
rho_p      = para(23,1);
rho_r      = para(24,1);
rho_nu     = para(25,1);
rho_z      = para(26,1);
rho_psi    = para(27,1);
rho_i      = para(28,1);
 rho_efp    = para(29,1);
 rho_nw     = para(30,1);

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
r_E    = exp(0.01*r_E_bar);

% define matrices of canonical system

neq  = 40;      %  Num of stable and unstable Variables
nex  = 10;       %  Num of Shock
nend = 10;       %  Num of Unstable Variables

GAM0j = zeros(neq,neq);
GAM1j = zeros(neq,neq);
C = zeros(neq,1);
PSI0j = zeros(neq,nex);
PPIj = zeros(neq,nend);

    GAM0j(1,1) = sigma/((theta/z_star-1.0)*((pii*theta)/r_n- ...
                 1.0))+(pii*sigma*theta^2)/(r_n*z_star*(theta/z_star-1.0...
                  )*((pii*theta)/r_n-1.0));
      GAM0j(1,2) = 1.0;
      GAM0j(1,18) = 1.0/((pii*theta)/r_n-1.0);
      GAM0j(1,28) = (sigma*theta)/(z_star*(theta/z_star-1.0)*( ...
                    (pii*theta)/r_n-1.0)) ;
      GAM0j(1,33) = -(pii*sigma*theta)/(r_n*(theta/z_star-1.0) ...
                  *((pii*theta)/r_n-1.0));
      GAM0j(1,39) = -(pii*theta)/(r_n*((pii*theta)/r_n-1.0));
      GAM0j(1,40) = -(pii*sigma*theta)/(r_n*(theta/z_star-1.0)...
                   *((pii*theta)/r_n-1.0));
      GAM0j(2,2) = 1.0;
      GAM0j(2,17) = -1.0;
      GAM0j(2,34) = -1.0;
      GAM0j(2,37) = 1.0;
      GAM0j(2,40) = sigma;
      GAM0j(3,2) = ((ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))/(ksi_w ...
                  *((kai*(lambda_w+1.0))/lambda_w+1.0)) ;
      GAM0j(3,3) = (pii*z_star)/r_n+((ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))...
                 /(ksi_w*((kai*(lambda_w+1.0))/lambda_w+1.0))+1.0;
      GAM0j(3,11) = (gamma_w*pii*z_star)/r_n+1.0;
      GAM0j(3,18) = -((ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))/(ksi_w... 
                   *((kai*(lambda_w+1.0))/lambda_w+1.0));
      GAM0j(3,20) = -1.0;
      GAM0j(3,28) = 1.0;
      GAM0j(3,29) = -(kai*(ksi_w-1.0)*((ksi_w*pii*z_star)/r_n-1.0))...
                    /(ksi_w*((kai*(lambda_w+1.0))/lambda_w+1.0));
      GAM0j(3,35) = -(pii*z_star)/r_n;
      GAM0j(3,37) = -(pii*z_star)/r_n;
      GAM0j(3,40) = -(pii*z_star)/r_n;
      GAM0j(4,4) = 1.0;
      GAM0j(4,6) = -(lambda_i+1.0)/(lambda_i-n_k+1.0);
      GAM0j(4,7) = (lambda_i+1.0)/(lambda_i-n_k+1.0)-1.0;
      GAM0j(4,16) = -(lambda_i+1.0)/(lambda_i-n_k+1.0);
      GAM0j(5,6) = (delta-1.0)/(psi*r_n);
      GAM0j(5,8) = -(delta-1.0)/(psi*r_n)-1.0;
      GAM0j(5,25) = 1.0;
      GAM0j(5,31) = 1.0;
      GAM0j(6,5) = 1.0;
      GAM0j(6,6) = -mu_E;
      GAM0j(6,7) = mu_E;
      GAM0j(6,16) = -mu_E;
      GAM0j(6,17) = -1.0;
      GAM0j(6,26) = -1.0;
      GAM0j(7,6) = ((delta-1.0)*(lambda_i+1.0))/(n_k*psi*r_E);
      GAM0j(7,7) = z_star/(eta*r_E);
      GAM0j(7,8) = -(((delta-1.0)/(psi*r_E)+1.0)*(lambda_i+1.0))/n_k;
      GAM0j(7,25) = (lambda_i+1.0)/n_k;
      GAM0j(7,27) = -1.0;
      GAM0j(7,28) = 1.0;
      GAM0j(8,3) = -1.0;
      GAM0j(8,8) = 1.0;
      GAM0j(8,9) = 1.0;
      GAM0j(8,25) = -1.0;
      GAM0j(8,28) = -1.0;
      GAM0j(8,29) = -1.0;
      GAM0j(9,6) = mu;
      GAM0j(9,8) = -mu;
      GAM0j(9,9) = 1.0;
      GAM0j(10,3) = alpha-1.0;
      GAM0j(10,8) = -alpha;
      GAM0j(10,10) = 1.0;
      GAM0j(11,10) = -((ksi_p-1.0)*((ksi_p*pii*z_star)/r_n-1.0))/ksi_p;
      GAM0j(11,11) = (gamma_p*pii*z_star)/r_n+1.0;
      GAM0j(11,21) = -1.0;
      GAM0j(11,37) = -(pii*z_star)/r_n;
      GAM0j(12,9) = -alpha*(phi_o_y+1.0);
      GAM0j(12,14) = 1.0;
      GAM0j(12,25) = alpha*(phi_o_y+1.0);
      GAM0j(12,28) = alpha*(phi_o_y+1.0);
      GAM0j(12,29) = (alpha-1.0)*(phi_o_y+1.0);
      GAM0j(13,1) = -c_y;
      GAM0j(13,14) = 1.0;
      GAM0j(13,15) = -i_y;
      GAM0j(13,19) = -1.0;
      GAM0j(14,9) = (delta+(psi*r_n)/pii-1.0)/(psi*z_star);
      GAM0j(14,15) = -(delta-1.0)/(psi*z_star)-1.0;
      GAM0j(14,16) = 1.0;
      GAM0j(14,22) = -(delta-1.0)/(psi*z_star)-1.0;
      GAM0j(14,25) = -(delta-1.0)/(psi*z_star);
      GAM0j(14,28) = -(delta-1.0)/(psi*z_star);
      GAM0j(15,6) = 1.0;
      GAM0j(15,15) = -inv_zeta-(inv_zeta*pii*z_star)/r_n;
      GAM0j(15,22) = 1.0;
      GAM0j(15,25) = -inv_zeta;
      GAM0j(15,28) = -inv_zeta;
      GAM0j(15,32) = -1.0;
      GAM0j(15,38) = (inv_zeta*pii*z_star)/r_n;
      GAM0j(15,40) = (inv_zeta*pii*z_star)/r_n;
      GAM0j(16,11) = phi_pi*(phi_r-1.0)*(1.0/4.0);
      GAM0j(16,12) = phi_pi*(phi_r-1.0)*(1.0/4.0);
      GAM0j(16,13) = phi_pi*(phi_r-1.0)*(1.0/4.0);
      GAM0j(16,14) = phi_y*(phi_r-1.0);
      GAM0j(16,17) = 1.0;
      GAM0j(16,23) = -1.0;
      GAM0j(17,5) = -1.0;
      GAM0j(17,30) = 1.0;
      GAM0j(17,37) = 1.0;
      GAM0j(18,1) = 1.0;
      GAM0j(19,2) = 1.0;
      GAM0j(20,3) = 1.0;
      GAM0j(21,31) = 1.0;
      GAM0j(22,8) = 1.0;
      GAM0j(23,11) = 1.0;
      GAM0j(24,15) = 1.0;
      GAM0j(25,18) = 1.0;
      GAM0j(26,28) = 1.0;
      GAM0j(27,18) = 1.0;
      GAM0j(28,19) = 1.0;
      GAM0j(29,20) = 1.0;
      GAM0j(30,21) = 1.0;
      GAM0j(31,22) = 1.0;
      GAM0j(32,23) = 1.0;
      GAM0j(33,24) = 1.0;
      GAM0j(34,25) = 1.0;
      GAM0j(35,32) = 1.0;
      GAM0j(36,26) = 1.0;
      GAM0j(37,27) = 1.0;
      GAM0j(38,24) = -1.0;
      GAM0j(38,25) = alpha/(alpha-1.0);
      GAM0j(38,28) = 1.0;
      GAM0j(39,12) = 1.0;
      GAM0j(40,13) = 1.0;

       GAM1j(1,1) = (sigma*theta)/(z_star*(theta/z_star-1.0)*...
                  (( pii*theta)/r_n-1.0));
      GAM1j(3,3) = 1.0;
      GAM1j(3,11) = gamma_w;
      GAM1j(5,6) = -1.0;
      GAM1j(7,6) = -(lambda_i+1.0)/n_k;
      GAM1j(7,7) = 1.0;
      GAM1j(8,16) = -1.0;
      GAM1j(11,11) = gamma_p;
      GAM1j(12,16) = alpha*(phi_o_y+1.0);
      GAM1j(14,16) = -(delta-1.0)/(psi*z_star);
      GAM1j(15,15) = -inv_zeta;
      GAM1j(16,13) = phi_pi*(phi_r-1.0)*(-1.0/4.0);
      GAM1j(16,14) = phi_y*(phi_r-1.0);
      GAM1j(16,17) = phi_r;
      GAM1j(18,33) = 1.0;
      GAM1j(19,34) = 1.0;
      GAM1j(20,35) = 1.0;
      GAM1j(21,30) = 1.0;
      GAM1j(22,36) = 1.0;
      GAM1j(23,37) = 1.0;
      GAM1j(24,38) = 1.0;
      GAM1j(25,39) = 1.0;
      GAM1j(26,40) = 1.0;
      GAM1j(27,18) = rho_b;
      GAM1j(28,19) = rho_g;
      GAM1j(29,20) = rho_w;
      GAM1j(30,21) = rho_p;
      GAM1j(31,22) = rho_nu;
      GAM1j(32,23) = rho_r;
      GAM1j(33,24) = rho_z;
      GAM1j(34,25) = rho_psi;
      GAM1j(35,32) = rho_i;
      GAM1j(36,26) = rho_efp;
      GAM1j(37,27) = rho_nw;
      GAM1j(39,11) = 1.0;
      GAM1j(40,12) = 1.0;

      PSI0j(27,1) = 1.0;
      PSI0j(28,2) = 1.0;
      PSI0j(29,3) = 1.0;
      PSI0j(30,4) = 1.0;
      PSI0j(31,5) = 1.0;
      PSI0j(32,6) = 1.0;
      PSI0j(33,7) = 1.0;
      PSI0j(34,8) = 1.0;
      PSI0j(35,9) = 1.0;
      PSI0j(36,10) = 1.0;
      PSI0j(37,11) = 1.0;

      PPIj(18,1) = 1.0;
      PPIj(19,2) = 1.0;
      PPIj(20,3) = 1.0;
      PPIj(21,4) = 1.0;
      PPIj(22,5) = 1.0;
      PPIj(23,6) = 1.0;
      PPIj(24,7) = 1.0;
      PPIj(25,8) = 1.0;
      PPIj(26,9) = 1.0;

% QZ(generalized Schur) decomposition by GENSYS
%[T1,TC,T0,TY,M,TZ,TETA,GEV,RC) = gensys(GAM0,GAM1,C,PSI0,PPI,1,1);
[T1,TC,T0,fmat,fwt,ywt,gev,RC,loose] = gensys(GAM0j,GAM1j,C,PSI0j,PPIj,1);
