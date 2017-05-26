function PlotPC_Mice_Model(data_mice, data_model, action, varargin)

color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black
BlockLabels = {'DA-L', 'DA-R','DA-Both', 'No DA'}; 
markerShape = {'o','s','v','^'};

% changing y-axis values
if any(data_mice(:,3)==-1)
  data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
end
if any(data_model(:,22)==-1)
  data_model(:,22) = (1 + data_model(:,22)) ./ 2;
end
if exist('action', 'var')
  action = (1+action)/2;
else
  action = data_model(:,22);
end

% Contrasts to evaluate P(R) at:
contrasts = unique(data_mice(:,2))';
if nargin >3
    if ~isempty(varargin{1})
      includeContrast = varargin{1};
    else
      includeContrast = ones(length(contrasts)*length(unique(data_mice(:,8))),1);        
    end
else
    includeContrast = ones(length(contrasts)*length(unique(data_mice(:,8))),1);
end
if nargin ==5
   plot_title = varargin{2}; 
else
    plot_title = 'Psychometric curves';
end

per_mice = nan(1,length(contrasts));
per_model = nan(1,length(contrasts));

niters = size(action,2);
b=0;

if nargin<6
  fig = figure('Position', [100, 100, 550, 450]);
  ax = axes('position',[0.09,0.1,0.86,0.88],'FontSize',16); hold on

else
  fig1 = varargin{3};
  fig = subplot(varargin{3}(1), varargin{3}(2), varargin{3}(3))
  ax = axes('position',[0.09,0.1,0.86,0.88],'FontSize',16); hold on
  
end
h = zeros(length(unique(data_mice(:,8))')+2,1);

for blockID = unique(data_mice(:,8))'
   b=b+1;
   c=1;
   for istim = unique(data_mice(:,2))'
      
      id=(data_mice(:,2)==istim & data_mice(:,8)==blockID);
      per_mice(1,c)= nanmean(data_mice(id,3)) ;
      
      per_model(c)=0;
      
      for iter = 1:niters
          id =(data_model(:,2)==istim & data_model(:,8)==blockID);
          per_model(c)= per_model(c) + nanmean(action( id,iter)) ;
      end
      per_model(c)=per_model(c)/niters;
      c=c+1;
   end
   plot_contrast = includeContrast((1:length(contrasts))+ (b-1)*length(contrasts) )==1;
  
   
   % PLOT
   % plot the mice raw data as points
  
   plot(contrasts(plot_contrast),per_mice(plot_contrast),...
      'color',color(blockID,:),...
      'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color(blockID,:),...
      'markerfacecolor',color(blockID,:),'linestyle','none',...
      'linewidth',1.2);
   plot(contrasts(plot_contrast),per_model(plot_contrast),'color',color(blockID,:),...
      'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color(blockID,:),...
      'linestyle','--','linewidth',1.2);

   h( b ) = plot(ax,NaN,NaN,'color',color(blockID,:),'marker',markerShape{blockID},'markersize',10,...
  'markerfacecolor',color(blockID,:),'markerfacecolor',color(blockID,:),'linestyle','none');
   ylabel(ax,'Fraction rightward choice')
   xlabel(ax,'Contrast')

end
   
% dummy plots for the legend

b=length(unique(data_mice(:,8))');
h(b+1) = plot(ax,NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
h(b+2) = plot(ax,NaN,NaN,'color','k','linestyle','--','linewidth',1.2);


legend(h, BlockLabels{ unique(data_mice(:,8))'},'Mouse','Model','Location','southeast');
title(plot_title)



end

