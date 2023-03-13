%/*******************************************************/
%/*                                                     */
%/* Computing Impulse responses based on                */
%/* parameters of DSGE models                           */
%/*                                                     */
%/*******************************************************/

% tic
% close all
% clear all
% clc

l = path;
% path('Impulse Responses',path);
path('Minnesota Prior',path);
path('Impulse Responses',path);
path('function',path);
addpath('DSGE_Prior');
addpath('gensys');
path('Data',path);
nvar = 6;
nshock=9;

[ZZ] = make_zz_NK();

nirf = 20; 
DSGE_nsim=100;
DSGE_nblock=10;
DSGE_nburnin_block= round(0.5*DSGE_nblock); 
runname = 'sw2011';
[para] = readpara(DSGE_nsim,DSGE_nblock,DSGE_nburnin_block,runname);

% 
 hirf = zeros(1,nvar*nshock*nirf);

 [ssirf, PP_NK ] = dsgeirf_NK(para,nirf,nvar,nshock,ZZ);
 hirf = reshape(ssirf,nvar*nshock*nirf,1)';
% end

% 
 irf_m_NK     = reshape(ssirf,nirf,nvar*nshock);


titlestr = {'Preference Shock: Z^b ' ,...
           'Exg-Goods Demand Shock: Z^g' ,...
		   'Wage Markup Shock: Z^w ' ,...
		   'Price Markup Shock: Z^p ' ,...
		   'Investment Price Shock: Z^{\nu}' ,...
		   'Monetary Policy Shock: Z^r' ,...
		   'Productivity Shock: Z^z ' ,... 
		   'Investment Shock: Z^{\psi}' ,...
           'Investment Price Shock: Z^{i}' ,...
		   'Financial Shock: Z^{efp}' ,...
		   'Net Worth Shock: Z^{nw}'};
%
ystr = {' -> Output' ,...
        ' -> Consumption' ,...
			' -> Investment' ,...  
			' -> Real Wage' ,...
			' -> Labor' ,...  
            ' -> Inflation' ,...
		    ' -> Investment Price' ,...
             ' -> Nominal Rate' ,...
             ' -> Real Borrowing', ...
             ' -> Loan Rate' };
             

% for i = 1:nshock
%    figure(3000+10*i)  
%   for j = 1:nvar
%     subplot(2, 3, j)
%     plot(1:nirf, irf_m(:,nvar*(i-1)+j),'b')
%     title(strcat(titlestr(i),ystr(j)))
%   end
% %   w = waitforbuttonpress;
% end

