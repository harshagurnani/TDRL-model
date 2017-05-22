function PlotCondPC_Mice_Model(data_mice, data_model, action, correct)

color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black

% changing y-axis values
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
data_model(:,22) = (1 + data_model(:,22)) ./ 2;
action = (1+action)/2;

xx = unique(data_mice(:,2))';
psycho_xx = min(xx):0.01:max(xx);
nn = nan(2,length(xx));
pp = nan(2,length(xx));

per_mice = nan(2,length(unique(data_mice(:,2))));
per_model = nan(2,length(unique(data_model(:,2))));




niters = size(action,2);
for blockID = unique(data_mice(:,8))'
   c=1;
   for istim = unique(data_mice(:,2))'
      
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          data_mice(1:end-1,3)==0 & data_mice(1:end-1,10)==1];
      per_mice(1,c)= nanmean(data_mice(id,3)) ;
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          data_mice(1:end-1,3)==1 & data_mice(1:end-1,10)==1];
      per_mice(2,c)= nanmean(data_mice(id,3)) ;
      
      
      per_model(1,c)=0;
      per_model(2,c)=0;
      
      for iter = 1:niters
          id =[false;data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==0 & correct(1:end-1,iter)==1];
          per_model(1,c)= per_model(1,c) + nanmean(action( id,iter)) ;
          id = [false; data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==1 & correct(1:end-1,iter)==1];
          per_model(2,c)= per_model(2,c) + nanmean(action( id,iter)) ;
      end
      per_model(1,c)=per_model(1,c)/niters;
      per_model(2,c)=per_model(2,c)/niters;
%      nn(1,c) = length(data_mice(data_mice(:,2)==istim & ...
%          data_mice(:,8)==blockID,2));
%       pp(1,c) = length(data_mice(data_mice(:,2)==istim & ...
%          data_mice(:,8)==blockID & data_mice(:,3)==1,2))/nn(1,c);
      
      c=c+1;
   end
   
%    psychoParams = mle_fit_psycho([xx ;nn ;pp],'erf_psycho');
%    
%    [~, V] = binostat(1,pp);  % variance of binomial distribution
   
   
   % PLOT
   % plot the mice raw data as points
   fig(blockID) = figure('Position', [100, 100, 550, 450]);
   ax(blockID) = axes(fig(blockID),'position',[0.09,0.1,0.86,0.88],'FontSize',16); hold on
   suptitle(sprintf('Psychomeric curve conditional to previous correct choice for Block %d', blockID))
   plot(ax(blockID),unique(data_mice(:,2))',per_mice(1,:),...
      'color',color(blockID,:),...
      'marker','o','markersize',10,'markeredgecolor','b',...
      'markerfacecolor','b','linestyle','none',...
      'linewidth',1.2);
   plot(ax(blockID),unique(data_mice(:,2))',per_mice(2,:),...
      'color',color(blockID,:),...
      'marker','o','markersize',10,'markeredgecolor','r',...
      'markerfacecolor','r','linestyle','none',...
      'linewidth',1.2);
   % plot plot error bars with sem across sessions
%    errorbar(ax,unique(data_mice(:,2))',per_mice,...
%       V ./ sqrt(length(unique(data_mice(:,9)))),'color',color(blockID,:),...
%       'marker','o','markersize',10,'markeredgecolor',color(blockID,:),...
%       'markerfacecolor',color(blockID,:),'linestyle','none',...
%       'linewidth',1.2)
%    
%    %plot the pscychometric function fitted to the mice data
%    plot(ax,psycho_xx, erf_psycho(psychoParams,psycho_xx),...
%       'color',color(blockID,:),'linestyle','-','linewidth',1.2);
   
   % plot the model
   plot(ax(blockID),unique(data_model(:,2))',per_model(1,:),'color','b',...
      'marker','o','markersize',10,'markeredgecolor','b',...
      'linestyle','--','linewidth',1.2);
   plot(ax(blockID),unique(data_model(:,2))',per_model(2,:),'color','r',...
      'marker','o','markersize',10,'markeredgecolor','r',...
      'linestyle','--','linewidth',1.2);

   ylabel(ax(blockID),'Fraction rightward choice')
   xlabel(ax(blockID),'Contrast')
   
   % dummy plots for the legend
   h = zeros(4,2);
   h(1,blockID) = plot(ax(blockID),NaN,NaN,'color','b','marker','o','markersize',10,...
      'markerfacecolor','b','linestyle','none');
   h(2,blockID) = plot(ax(blockID),NaN,NaN,'color','r','marker','o','markersize',10,...
      'markerfacecolor','r','linestyle','none');
   h(3,blockID) = plot(ax(blockID),NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
   h(4,blockID) = plot(ax(blockID),NaN,NaN,'color','k','linestyle','--','linewidth',1.2);
   
   
   legend(h(:,blockID), 'After correct L','After correct R','Mouse','Model','Location','southeast');
   
   
end



end

