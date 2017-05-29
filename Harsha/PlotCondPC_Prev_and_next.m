function PlotCondPC_Prev_and_next(data_mice, varargin)

color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black
blockLabels ={'DA-L', 'DA-R', 'DA-Both', 'No DA'};
markerShapes={'o', 's', 'v', '^'};

% changing y-axis values
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;

contrasts = unique(data_mice(:,2))';
per_mice_next = nan(2,length(contrasts));
per_mice_prev = nan(2,length(contrasts));

if nargin >1
    includeContrast = varargin{1}; 
else
    includeContrast = ones(length(contrasts)*length(unique(data_mice(:,8))),1);
end

if nargin>2
    blocks = varargin{2};
else
    blocks=unique(data_mice(:,8))';
end

fig = figure;%varargin{3};
ax = axes(fig,'position',[0.09,0.1,0.86,0.88],'FontSize',16); hold on 


b=0;
for blockID = blocks
   b=b+1;
   c=1;
   for istim = unique(data_mice(:,2))'
      
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          data_mice(1:end-1,3)==0 & data_mice(1:end-1,10)==1];
      per_mice_next(1,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & ...
          data_mice(2:end,3)==0 & data_mice(2:end,10)==1; false];
      per_mice_prev(1,c)= nanmean(data_mice(id,3)) ;
      
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          data_mice(1:end-1,3)==1 & data_mice(1:end-1,10)==1];
      per_mice_next(2,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & ...
          data_mice(2:end,3)==1 & data_mice(2:end,10)==1; false];
      per_mice_prev(2,c)= nanmean(data_mice(id,3)) ;
     
      c=c+1;
   end

   plot_contrast = includeContrast((1:length(contrasts))+ (b-1)*length(contrasts) )==1;
   
   %% PLOT
   % plot the mice raw data as points
%    if nargin<=6
%    fig(blockID) = figure('Position', [100, 100, 550, 450]);
%    ax = axes(fig(blockID),'position',[0.09,0.1,0.86,0.88],'FontSize',16); hold on   
%    suptitle(sprintf('Psychomeric curve conditional to previous correct choice for Block %s',blockLabels{blockID})) 
%    end

   %plot
   plot(ax,contrasts(plot_contrast),per_mice_next(1,plot_contrast),...
      'color','b',...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor','b',...
      'markerfacecolor','b','linestyle','-',...
      'linewidth',1.2);
   plot(ax,contrasts(plot_contrast),per_mice_next(2,plot_contrast),...
      'color','r',...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor','r',...
      'markerfacecolor','r','linestyle','-',...
      'linewidth',1.2);

   
   % plot the model
   plot(ax,contrasts(plot_contrast),per_mice_prev(1,plot_contrast),'color','b',...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor','b',...
      'linestyle','--','linewidth',1.2);
   plot(ax,contrasts(plot_contrast),per_mice_prev(2,plot_contrast),'color','r',...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor','r',...
      'linestyle','--','linewidth',1.2);

   ylabel(ax,'Fraction rightward choice')
   xlabel(ax,'Contrast')
   
   % dummy plots for the legend
   
%    if nargin <=6
%    h = zeros(4,1);
%    h(1) = plot(ax,NaN,NaN,'color','b','marker',markerShapes{blockID},'markersize',10,...
%       'markerfacecolor','b','linestyle','none');
%    h(2) = plot(ax,NaN,NaN,'color','r','marker',markerShapes{blockID},'markersize',10,...
%       'markerfacecolor','r','linestyle','none');
%    h(3) = plot(ax,NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
%    h(4) = plot(ax,NaN,NaN,'color','k','linestyle','--','linewidth',1.2);
%    
%    
%    legend(h, 'After correct L','After correct R','Mouse','Model','Location','southeast');
%    end
   
end
%    
% if nargin>6
   h = zeros(4+length(blocks),1);
   h(1) = plot(ax,NaN,NaN,'color','b','marker','o','markersize',10,...
      'markerfacecolor','b','linestyle','none');
   h(2) = plot(ax,NaN,NaN,'color','r','marker','o','markersize',10,...
      'markerfacecolor','r','linestyle','none');
   h(3) = plot(ax,NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
   h(4) = plot(ax,NaN,NaN,'color','k','linestyle','--','linewidth',1.2);
   b=1;
   for blockID = blocks
      h(4+b) =  plot(ax,NaN,NaN,'color','k','marker',markerShapes{blockID},'markersize',10,...
      'markerfacecolor','k','linestyle','none');
      b=b+1;
   end
   legend(h, 'Condtioned to correct L','Condtioned to correct R','After Correct choice','Before correct choice',blockLabels{blocks},'Location','southeast');
% end


end

