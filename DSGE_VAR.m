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

% tic
% close all
% clear all
% clc

l = path;
% path('Impulse Responses',path);
path('Minnesota Prior',path);
path('function',path);
addpath('DSGE_Prior');
addpath('gensys');
path('Data',path);


%=========================================================================
%         GENERATE DUMMY OBSERVATIONS FROM MINNESOTA PRIOR 
%=========================================================================

%   lambda = 1;
%   mu = 1;

% p          = 2;                 % Number of lags in the VAR

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
% Sigma     = 1/(1+lambda)*((lambda*YYdum+YY)-(lambda*XYdum+XY)'...
%             * inv(lambda*XXdum+XX)*(lambda*XYdum+XY) );
        
Sigma     = 1/(1+lambda)*(YY-(XY)'*inv(XX)*(XY) );         
Sigma    = 0.5*(Sigma + Sigma')  ;      

% Matrices for collecting draws from Posterior Density

Sigmap    = zeros(nsim,n,n);
Phip      = zeros(nsim,n*p+1,n);
largeeig  = zeros(nsim,1);
counter   = 0;

%=========================================================================
%            DRAWS FROM POSTERIOR DENSITY (DIRECT SAMPLING)
%=========================================================================
% disp('                                                                  ');
% disp('        BAYESIAN ESTIMATION OF VAR: DIRECT SAMPLING...            ');
% disp('                                                                  ');

for j=1:nsim

    
    % Draws from the density Sigma | Y
    %  逆ウィシャート分布による分散・共分散分布のサンプルの生成
    sigma   = iwishrnd(T*Sigma,T*(1+lambda)-n*p-1);
    
    % Draws from the density vec(Phi) |Sigma(j), Y
    %  多変量正規分布による係数のサンプルの生成    
    
    XX1 = 0.5*((XX)'+(XX));    
    phi_new = mvnrnd(reshape(Phi_tilde,n*(n*p+1),1),kron(sigma,inv(T*(XX1))));
    
    % Rearrange vec(Phi) into Phi
    
    Phi     = reshape(phi_new,n*p+1,n);
        
    Sigmap(j,:,:) = sigma;
    Phip(j,:,:)   = Phi;
    
    Phi = Phi(1:n*p,:);

    
    % Calculating largest eigenvalue of Companion form
    
     F(1:n,1:n*p)    = Phi';

     eigen           = eig(F);
     eigen           = max(eigen);
     largeeig(j)     = abs(eigen);
     counter         = counter +1; 
     
    if mod(j,1000)==0
%        disp(['         DRAW NUMBER:   ', num2str(j)]);
%        disp('                                                                  ');
%        disp(['     REMAINING DRAWS:   ', num2str(nsim-j)]);
       disp('                                                                  ');
    end
     
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



% %=========================================================================
% %           FIGURE 1: LARGEST EIGENVALUE (Companion Form)
% %=========================================================================
% 
% pnames = strvcat( 'Largest Eigenvalue (Recursive Average)',...
%     'Largest Eigenvalue (Posterior Marginal Density)');
% 
% figure('Position',[20,20,900,600],'Name',...
%     'Largest Eigenvalue (Companion Form)','Color','w')
% 
% rmean = zeros(nsim,1);
% 
% for i=1:nsim
%     rmean(i) = mean(largeeig(1:i));
% end
% 
% subplot(1,2,1), plot(rmean,'LineStyle','-','Color','b',...
%         'LineWidth',2.5), hold on
% title(pnames(1,:),'FontSize',13,'FontWeight','bold');
% 
% [density,x]  = ksdensity(largeeig(nburn:end));
% 
% subplot(1,2,2), plot(x,density,'LineStyle','-','Color','b',...
%         'LineWidth',2.5), hold on
% 
% title(pnames(2,:),'FontSize',13,'FontWeight','bold');
% 
% path=l;

% disp(['         ELAPSED TIME:   ', num2str(toc)]);

% elapsedtime=toc;
