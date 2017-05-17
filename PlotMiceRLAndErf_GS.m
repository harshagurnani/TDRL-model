function PlotMiceRLAndErf_GS(data_mice, data_model)

color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black

% changing y-axis values
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
data_model(:,16) = (1 + data_model(:,16)) ./ 2;


xx = unique(data_mice(:,2))';
psycho_xx = min(xx):0.01:max(xx);
nn = nan(1,length(xx));
pp = nan(1,length(xx));

per_mice = nan(1,length(unique(data_mice(:,2))));
per_model = nan(1,length(unique(data_model(:,2))));


fig = figure('Position', [100, 100, 550, 450]);
ax = axes(fig,'position',[0.09,0.1,0.86,0.88],'FontSize',16); hold on

for blockID = unique(data_mice(:,8))'
   c=1;
   for istim = unique(data_mice(:,2))'
      
      per_mice(c)= nanmean(data_mice(data_mice(:,2)==istim & ...
         data_mice(:,8)==blockID,3)) ;
      per_model(c)= nanmean(data_model(data_model(:,2)==istim & ...
         data_model(:,8)==blockID,16)) ;
      
      nn(c) = length(data_mice(data_mice(:,2)==istim & ...
         data_mice(:,8)==blockID,2));
      pp(c) = length(data_mice(data_mice(:,2)==istim & ...
         data_mice(:,8)==blockID & data_mice(:,3)==1,2))/nn(c);
      
      c=c+1;
   end
   
   psychoParams = mle_fit_psycho([xx ;nn ;pp],'erf_psycho');
   
   [~, V] = binostat(1,pp);  % variance of binomial distribution
   
   
   % PLOT
   % plot the mice raw data as points
   plot(ax,unique(data_mice(:,2))',per_mice,...
      'color',color(blockID,:),...
      'marker','o','markersize',10,'markeredgecolor',color(blockID,:),...
      'markerfacecolor',color(blockID,:),'linestyle','none',...
      'linewidth',1.2);
   
   % plot plot error bars with sem across sessions
   errorbar(ax,unique(data_mice(:,2))',per_mice,...
      V ./ sqrt(length(unique(data_mice(:,9)))),'color',color(blockID,:),...
      'marker','o','markersize',10,'markeredgecolor',color(blockID,:),...
      'markerfacecolor',color(blockID,:),'linestyle','none',...
      'linewidth',1.2)
   
   %plot the pscychometric function fitted to the mice data
   plot(ax,psycho_xx, erf_psycho(psychoParams,psycho_xx),...
      'color',color(blockID,:),'linestyle','-','linewidth',1.2);
   
   % plot the model
   plot(ax,unique(data_model(:,2))',per_model,'color',color(blockID,:),...
      'marker','o','markersize',10,'markeredgecolor',color(blockID,:),...
      'linestyle','--','linewidth',1.2);
   
   
end

   ylabel(ax,'Fraction rightward choice')
   xlabel(ax,'Contrast')
   
   % dummy plots for the legend
   h = zeros(4,1);
   h(1) = plot(ax,NaN,NaN,'color','b','marker','o','markersize',10,...
      'markerfacecolor','b','linestyle','none');
   h(2) = plot(ax,NaN,NaN,'color','r','marker','o','markersize',10,...
      'markerfacecolor','r','linestyle','none');
   h(3) = plot(ax,NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
   h(4) = plot(ax,NaN,NaN,'color','k','linestyle','--','linewidth',1.2);
   
   
   legend(h, 'Left DA','Right DA','Mouse','Model','Location','southeast');


end

