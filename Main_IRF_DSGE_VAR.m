%% DOCUMENT TITLE
% INTRODUCTORY TEXT
%%
%==========================================================================
%                       VAR WITH MINNESOTA PRIOR                      
%
%                      See README_VAR for details
%
%
%==========================================================================


%=========================================================================
%                             HOUSEKEEPING
%=========================================================================

tic
close all
clear all
clc

l = path;
% path('Impulse Responses',path);
path('Minnesota Prior',path);
path('Impulse Responses',path);
path('function',path);
addpath('DSGE_Prior');
addpath('gensys');
path('Data',path);


%=========================================================================
%         GENERATE DUMMY OBSERVATIONS FROM MINNESOTA PRIOR 
%=========================================================================

H    = 12;    % Number of Periods for the computation of IRFs

shock = 6  % monetary 6:

   lambda = 0.4 %0.1; %4;
    mu = 0.5 ;

p  = 2;                 % Number of lags in the VAR

% ダミー変数型の事前分布のハイパーパラメータの設定 
% vm_spec

% 事前分布に基づくダミー変数の作成
DSGE_dummy

%=========================================================================
%     DEFINITION OF DATA, LAG STRUCTURE AND POSTERIOR SIMULATION
%=========================================================================

% [Tdummy,n] = size(YYdum);
[Tobs,n]   = size(YYact);
X          = [XXact];
Y          = [YYact];
n          = n;                 % Number of variables in the VAR
p          = p;                 % Number of lags in the VAR
T          = Tobs;


XX = lambda*XXdum + (X'*X)/T;
XY = lambda*XYdum + (X'*Y)/T;
YY = lambda*YYdum + (Y'*Y)/T;

% サンプル数
nsim       = 1000;             % Number of draws from Posterior Density
% バーンインのサンプル数
nburn      = 0.2*nsim;          % Number of draws to discart

F          = zeros(n*p,n*p);    % Matrix for Companion Form
I          = eye(n);

for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end


%=========================================================================
%               OLS ESTIMATOR FOR PHI AND SSR (SIGMA)
%=========================================================================

% 係数の推定値
% Phi_tilde = inv(lambda*XXdum+XX)*(lambda*XYdum+XY);
Phi_tilde = inv(XX)*(XY);

% 分散共分散行列の推定値    
Sigma     = 1/(1+lambda)*(YY-(XY)'*inv(XX)*(XY) );           
Sigma    = 0.5*(Sigma + Sigma')  ;      

% Matrices for collecting draws from Posterior Density

Sigmap    = zeros(nsim,n,n);
Phip      = zeros(nsim,n*p+1,n);

IRF_Pref  = zeros(nsim,n,H+1);
IRF_MP    = zeros(nsim,n,H+1);

largeeig  = zeros(nsim,1);
counter   = 0;

DSGE_irf
DSGE_irf_NK

DSGE_nsim=100;
DSGE_nblock=10;
DSGE_nburnin_block= round(0.5*DSGE_nblock); 
nuse = 1;
runname = 'sw2011';
iparasim = strcat(runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);
indseq0 = mod(1:(DSGE_nblock-DSGE_nburnin_block)*DSGE_nsim, nuse);
indseq = (indseq0 ~= 0);
parasim_post_NK  = delif(fhpara,indseq);

runname = 'FF2011';
iparasim = strcat(runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);
indseq0 = mod(1:(DSGE_nblock-DSGE_nburnin_block)*DSGE_nsim, nuse);
indseq = (indseq0 ~= 0);
parasim_post_FF  = delif(fhpara,indseq);


%=========================================================================
%            DRAWS FROM POSTERIOR DENSITY (DIRECT SAMPLING)
%=========================================================================
% disp('                                                                  ');
% disp('        BAYESIAN ESTIMATION OF VAR: DIRECT SAMPLING...            ');
% disp('                                                                  ');

for j=1:nsim
    
%     dsge_type = 1; % 'NK_model';
%     para_NK = (parasim_post_NK(j,1:end))' ;
%        [XXdum1, XYdum1, YYdum1] = dsge_prior(dsge_type,para_NK,nlags_);
%   
%     dsge_type = 2; % 'FF_model';
%     para_FF = (parasim_post_FF(j,1:end))' ;
%        [XXdum2, XYdum2, YYdum2] = dsge_prior(dsge_type,para_FF,nlags_);
%     %%
%     XXdum = mu*XXdum1+ (1-mu)*XXdum2;
%     XYdum = mu*XYdum1+ (1-mu)*XYdum2;
%     YYdum = mu*YYdum1+ (1-mu)*YYdum2;
%     XX = lambda*XXdum + (X'*X)/T;
%     XY = lambda*XYdum + (X'*Y)/T;
%     YY = lambda*YYdum + (Y'*Y)/T;
%     %%
%     Phi_tilde = inv(XX)*(XY);
%     Sigma     = 1/(1+lambda)*((YY)-(XY)'*inv(XX)*(XY) );           
%     Sigma    = 0.5*(Sigma + Sigma')  ;      

    
    % Draws from the density Sigma | Y
    %  逆ウィシャート分布による分散・共分散分布のサンプルの生成
    sigma   = iwishrnd(T*Sigma,T*(1+lambda)-n*p-1);
    
    % Draws from the density vec(Phi) |Sigma(j), Y
    %  多変量正規分布による係数のサンプルの生成    
    
    XX1 = 0.5*((XX)'+(XX));    
    phi_new = mvnrnd(reshape(Phi_tilde,n*(n*p+1),1),kron(sigma,inv(T*(XX1)))); 
   
    % Rearrange vec(Phi) into Phi
    
    Phi     = reshape(phi_new,n*p+1,n);
    
    % Computation of Impulse Responce using sign restrictions for the first
    % two periods; 
    
    shock =1;
    IR_Pref          = irf_DSGE_VAR(Phi,sigma,H+1,2,PP_FF,PP_NK,mu,shock);
    shock =6;
    IR_MP            = irf_DSGE_VAR(Phi,sigma,H+1,2,PP_FF,PP_NK,mu,shock);
    
        
    Sigmap(j,:,:) = sigma;
    Phip(j,:,:)   = Phi;
    
    Phi = Phi(1:n*p,:);
    IRF_Pref(j,:,:)  = IR_Pref;
    IRF_MP(j,:,:)    = IR_MP;
    
    % Calculating largest eigenvalue of Companion form
    
     F(1:n,1:n*p)    = Phi';

     eigen           = eig(F);
     eigen           = max(eigen);
     largeeig(j)     = abs(eigen);
     counter         = counter +1; 
     
%     if mod(j,100)==0
%        disp(['         DRAW NUMBER:   ', num2str(j)]);
%        disp('                                                                  ');
%        disp(['     REMAINING DRAWS:   ', num2str(nsim-j)]);
%        disp('                                                                  ');
%     end
     
end

%=========================================================================
%             Posterior Distributions of Parameters
%=========================================================================

%   Cal_Post;


%=========================================================================
%                        MARGINAL DATA DENSITY
%=========================================================================

 vm_mdd
% 
 mdd = real(lnpYY);               % Marginal Data Density
% 
 disp(['Marginal Likelihood is', num2str(mdd) ]);



%=========================================================================
% Compute mean and 90% pointwise credible sets (HPD) for Impulse Responces
%=========================================================================

[irmean_MP,ir05_MP,ir95_MP] = impulsemoment(IRF_MP);
[irmean,ir05,ir95] = impulsemoment(IRF_Pref);


%=========================================================================
%           Figure 1: Impulse Responce Functions (Posterior means and 90%
%                     Credible Sets)
%=========================================================================

pnames = strvcat( 'Output','Cons','inv','real wage','Inflation', 'Interest Rate');

figure('Position',[20,20,900,600],'Name',...
        'IRFs to a Preference Shock','Color','w')
j= 1; i=2;
scale = 10;
scale_d = scale*irmean(i,1)/irf_m(1,nvar*(j-1)+i);
for i=1:6
subplot(3,2,i),  
hold on
   h=area([scale*ir05(i,:)' scale*(ir95(i,:)-ir05(i,:))' ] );
      set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
      set(h(2),'FaceColor',[0.65 1 1])      
      set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value      
     
plot(scale*irmean(i,:),'LineStyle',':','Color','k','LineWidth',2.5),  

plot(1:H+1, scale_d*irf_m(1:H+1,nvar*(j-1)+i),'LineStyle','-','Color','r',...
        'LineWidth',2.5),       
plot(1:H+1, scale_d*irf_m_NK(1:H+1,nvar*(j-1)+i),'LineStyle','--','Color','b',...
        'LineWidth',2.5),        
hold off
title(pnames(i,:),'FontSize',13,'FontWeight','bold');
end
legend('','95% Band','MDSGE-VAR','DSGE(FF)','DSGE(NK)')


%%  Monetary Policy Scock
figure('Position',[20,20,900,600],'Name',...
        'IRFs to a Monetary Policy Shock','Color','w')
j= 6; % shock
scale = 10;
scale_d = scale*irmean_MP(j,1)/irf_m(1,nvar*(j-1)+j);
for i=1:6
subplot(3,2,i),  
hold on
   h=area([scale*ir05_MP(i,:)' scale*(ir95_MP(i,:)-ir05_MP(i,:))' ] );
      set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
      set(h(2),'FaceColor',[0.65 1 1])      
      set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value      
     
plot(scale*irmean_MP(i,:),'LineStyle',':','Color','k','LineWidth',2.5),  

plot(1:H+1, scale_d*irf_m(1:H+1,nvar*(j-1)+i),'LineStyle','-','Color','r',...
        'LineWidth',2.5),       
plot(1:H+1, scale_d*irf_m_NK(1:H+1,nvar*(j-1)+i),'LineStyle','--','Color','b',...
        'LineWidth',2.5),        
hold off
title(pnames(i,:),'FontSize',13,'FontWeight','bold');
end
legend('','95% Band','MDSGE-VAR','DSGE(FF)','DSGE(NK)')


path = l;

disp(['         ELAPSED TIME:   ', num2str(toc)]);

elapsedtime=toc;


% disp(['         ELAPSED TIME:   ', num2str(toc)]);

% elapsedtime=toc;
