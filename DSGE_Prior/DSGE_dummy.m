%******************************************************************/
%* filename :    vm_dummy.m                                       */
%* Description:  Generates dummy observations                     */
%*               for a Minnesota Prior                            */
%******************************************************************/

% close all; 
% clc;

nlags_   = p;   
T0       = 4;               % size of pre-sample */

%******************************************************** 
% Import data series                                    *
%*******************************************************/
vm_loaddata

nv      = size(YY,2);     %* number of variables */
nobs    = size(YY,1)-T0;  %* number of observations */

%******************************************************** 
% Dummy Observations                                    *
%*******************************************************/

%** Obtain mean and standard deviation from expandend pre-sample data
 
YY0     =   YY(1:T0+16,:);  
ybar    =   mean(YY0)';      
sbar    =   std(YY0)'; 
premom  =   [ybar sbar];

% Generate matrices with dummy observations
% hyp = [tau; d; w; lambda; mu];

% [YYdum, XXdum, breakss] = varprior_h(nv,nlags_,nex_,hyp,premom);
DSGE_nsim=100;
DSGE_nblock=10;
DSGE_nburnin_block= round(0.5*DSGE_nblock); 
runname = 'sw2011';
[para] = readpara(DSGE_nsim,DSGE_nblock,DSGE_nburnin_block,runname);

dsge_type = 1; % 'NK_model';
 [XXdum1, XYdum1, YYdum1] = dsge_prior(dsge_type,para,nlags_);
 
DSGE_nsim=100;
DSGE_nblock=10;
DSGE_nburnin_block= round(0.5*DSGE_nblock); 
runname = 'FF2011';
[para] = readpara(DSGE_nsim,DSGE_nblock,DSGE_nburnin_block,runname); 
 
dsge_type = 2; % 'FF_model';
 [XXdum2, XYdum2, YYdum2] = dsge_prior(dsge_type,para,nlags_);
 
 XXdum = mu*XXdum1+ (1-mu)*XXdum2;
 XYdum = mu*XYdum1+ (1-mu)*XYdum2;
 YYdum = mu*YYdum1+ (1-mu)*YYdum2;
 
% Actual observations

YYact = YY(T0+1:T0+nobs,:);
ti = ti(T0+1:T0+nobs,:);
XXact = zeros(nobs,nv*nlags_);
i = 1;

while (i <= nlags_)
    XXact(:,(i-1)*nv+1:i*nv) = YY(T0-(i-1):T0+nobs-i,:);
    i = i+1;
end

% last column of XXact = constant
XXact = [XXact ones(nobs,1)];

% store observations
% save data YYdum XXdum YYact XXact
