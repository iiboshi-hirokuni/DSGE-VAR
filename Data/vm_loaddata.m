%Data set:
%Output, Inflation, Interest Rates, InvVelocity  

%---------------------------------------------------------------------
% loading DATA
%---------------------------------------------------------------------

% datapath = './data/';
datafilename = 'kaihatsu_kurozumi_data6.csv';

% series_YT = csvread(strcat(datapath, datafilename), 1, 1);
series_YT = csvread(strcat(datafilename), 1, 1);

nobs = size(series_YT, 1);  % number of observations 
yy_m = mean(series_YT, 1);
YY = series_YT;
ti=linspace(1981,0.2, size(YY,1))';


% load datampshock.txt
% YY=datampshock;
% ti=linspace(1964,0.2, size(YY,1))';
% nobs=size(YY,1);

%Convert velocity into real money balances
% YY(:,4)=YY(:,4)+YY(:,1);
