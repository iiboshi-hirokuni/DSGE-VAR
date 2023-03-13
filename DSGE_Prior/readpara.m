
function [para] = readpara(nsim,nblock,nburnin_block,runname);

nuse    = 10; % use every nuse observation
nstate = 40;

% resupath ='';
% runname = 'FF2011';
iparasim = strcat(runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);

indseq0 = mod(1:(nblock-nburnin_block)*nsim, nuse);
indseq = (indseq0 ~= 0);
parasim_post  = delif(fhpara,indseq);
nasim = size(parasim_post,1);

para = mean(parasim_post(:,1:end))' ;