% Calculating of Posterior estimates 

a  =0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
rate = (1-a)/2;


% calculation of posterior estimates of parameters
%

sort_para=zeros( nsim-nburn, n*p+1, n );

for k = 1:1:n 
  for i=1:1:n*p+1
     sort_para(:,i,k) = sort(Phip(nburn+1:end,i,k),1);
  end
end

 for l=1:nv
%  The Fisrt equation
% j=1;
  para_low = sort_para(round((nsim-nburn)*rate),:,l); 
  para_up  = sort_para(round((nsim-nburn)*(1-rate)),:,l);

% calculating inefficiency
para_save = Phip(nburn:end,:,l);
cal_inefficiency

disp( [ num2str(l) '-th equation'] )
disp('    mean     [ lower    Upper ] inefficiency');
  for i=1:n*p+1
     disp([ mean(Phip(nburn:end,i,l)) para_low(i)  para_up(i) inefficiency(i) ] );
  end
 end

% for l=1:nv
% figure(1000+100*l)
% %     'Posteror Estimates of Coeffcients in GDP eq','Color','w')
% for i =1:n*p+1
%  subplot(3,3,i),    
%      [density,x1]  = ksdensity(Phip(nburn:end,i,l));
%     plot(x1,density,'LineStyle','-','Color','b',...
%         'LineWidth',2.5);
%      title( ['phi'  num2str(l)  num2str(i) ],'FontSize',12 );
% end     
% end

% 
% figure('Name',...
%     'Posteror Estimates of Covariance Matrix','Color','w')
% %figure(5);
% for j =1:n
%  for i =1:n
%  subplot(6,6,i+6*(j-1)),    
%      [density,x1]  = ksdensity(Sigmap(nburn:end,i,j));
%     plot(x1,density,'LineStyle','-','Color','b',...
%         'LineWidth',2.5);
%      title( 'Sigma','FontSize',12 );
%  end  
% end




