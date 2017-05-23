function PlotQValues( data_model, action, correct, QL, QR)

color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black

% changing y-axis values
data_model(:,22) = (1 + data_model(:,22)) ./ 2;
action = (1+action)/2;

per_model = nan(5,length(unique(data_model(:,2))));




niters = size(action,2);
for blockID = unique(data_model(:,8))'
   c=1;
   for istim = unique(data_model(:,2))'
    
      
      per_model(:,c)=0;
      
      for iter = 1:niters
          % QR-QL for diff stimuli
          id =[data_model(:,2)==istim & data_model(:,8)==blockID];
          per_model(1,c)= per_model(1,c) + nanmean(QR( id,iter) - QL( id,iter)) ;
          
          % QR-QL for diff stimuli after correct L
          id =[false;data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==0 & correct(1:end-1,iter)==1];
          per_model(2,c)= per_model(2,c) + nanmean(QR( id,iter) - QL( id,iter)) ;
          % QR-QL for diff stimuli after correct R
          id = [false; data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==1 & correct(1:end-1,iter)==1];
          per_model(3,c)= per_model(3,c) + nanmean(QR( id,iter) - QL( id,iter)) ;
          
          % QR-QL for diff stimuli after incorrect L
          id =[false;data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==0 & correct(1:end-1,iter)==0];
          per_model(4,c)= per_model(4,c) + nanmean(QR( id,iter) - QL( id,iter)) ;
          % QR-QL for diff stimuli after incorrect R
          id = [false; data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==1 & correct(1:end-1,iter)==0];
          per_model(5,c)= per_model(5,c) + nanmean(QR( id,iter) - QL( id,iter)) ;
      end
      per_model(:,c)=per_model(:,c)/niters;
      
      c=c+1;
   end

   %% PLOT
   % plot the mice raw data as points
   fig(blockID) = figure('Position', [100, 100, 550, 450]);
   ax(blockID) = axes(fig(blockID),'position',[0.09,0.1,0.86,0.88],'FontSize',16); hold on
   suptitle(sprintf('Psychomeric curve conditional to previous correct choice for Block %d', blockID))
   plot(ax(blockID),unique(data_model(:,2))',per_model(1,:),...
      'color','k',...
      'marker','o','markersize',10,'markeredgecolor','k',...
      'markerfacecolor','k','linestyle','none',...
      'linewidth',2);
   plot(ax(blockID),unique(data_model(:,2))',per_model(2,:),...
      'color','b',...
      'marker','o','markersize',10,'markeredgecolor','b',...
      'markerfacecolor','b','linestyle','none',...
      'linewidth',1.2);
   plot(ax(blockID),unique(data_model(:,2))',per_model(3,:),'color','r',...
      'marker','o','markersize',10,'markeredgecolor','r','markerfacecolor','r',...
      'linestyle','none','linewidth',1.2);
   plot(ax(blockID),unique(data_model(:,2))',per_model(4,:),'color','b',...
      'marker','o','markersize',10,'markeredgecolor','b',...
      'linestyle','--','linewidth',1.2);

   plot(ax(blockID),unique(data_model(:,2))',per_model(5,:),'color','r',...
      'marker','o','markersize',10,'markeredgecolor','r',...
      'linestyle','--','linewidth',1.2);
   ylabel(ax(blockID),'Fraction rightward choice')
   xlabel(ax(blockID),'Contrast')
   
   % dummy plots for the legend
   h = zeros(5,2);
   h(1,blockID) = plot(ax(blockID),NaN,NaN,'color','k','linestyle','-','linewidth',2,'marker','o','markerfacecolor','k','markersize',10,'markeredgecolor','k');
   h(2,blockID) = plot(ax(blockID),NaN,NaN,'color','b','marker','o','markersize',10,...
      'markerfacecolor','b','linestyle','none');
   h(3,blockID) = plot(ax(blockID),NaN,NaN,'color','r','marker','o','markersize',10,...
      'markerfacecolor','r','linestyle','none');
   h(4,blockID) = plot(ax(blockID),NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
   h(5,blockID) = plot(ax(blockID),NaN,NaN,'color','k','linestyle','--','linewidth',1.2);
   
   
   legend(h(:,blockID), 'All trials in block','After Left','After Right','After Correct','After Incorrect','Location','southeast');
   
   
end



end

