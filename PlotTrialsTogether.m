function PlotTrialsTogether(Data,data_model)
% nice plot

data_mice = Data.data;

% changing y-axis values
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
data_model(:,22) = (1 + data_model(:,22)) ./ 2;


plotLen = length(data_mice);
% if length(data_mice)<plotLen
%    
%    plotLen = length(data_mice);
%    
% end
plotLenVec = cat(1,ones(plotLen,1),zeros(length(data_mice)-plotLen,1));

mouse = data_mice(plotLenVec==1,3) ;
model = data_model(plotLenVec==1,22);

% smoothmouse = mouse*0.96 + 0.02;
% smoothmodel = model*0.96 + 0.02;

smoothmouse = smooth(mouse,7);
smoothmodel = smooth(model,7);


block1 = data_mice(plotLenVec==1,8);
block1(data_mice(:,8)==1) = -0.03;
block1(data_mice(:,8)~=1) = NaN;

block2 = data_mice(plotLenVec==1,8);
block2(data_mice(:,8)==2) = 1.03;
block2(data_mice(:,8)~=2) = NaN;
%block2 = data_mice(plotLenVec==1 & data_mice(:,8)==2,8) - 1;
% block1 = -0.03;
%block2 = 1.03;

% blockID = data_mice(plotLenVec==1,8) - 1;
% blockID(blockID==1) = 1.03;
% blockID(blockID==0) = -0.03;


fig = figure('Position', [100, 100, 1200, 350]);
ax = axes(fig,'position',[0.05,0.14,0.93,0.86],'fontsize',16); hold on

% plot(ax,mouse,'marker','o','MarkerEdgeColor','r','markerfacecolor',...
%    'r','markersize',1.5,'linestyle','none');

plot(ax,block1,'marker','o','MarkerEdgeColor','b','markerfacecolor',...
   'b','markersize',2,'linestyle','none');

plot(ax,block2,'marker','o','MarkerEdgeColor','r','markerfacecolor',...
   'r','markersize',2,'linestyle','none');

h(1) = plot(ax,smoothmouse,'linestyle','-','linewidth',1.2);

h(2) = plot(ax,smoothmodel,'color',[0 0 0],'linestyle','-','linewidth',1.2);

h(3) = plot(ax,NaN,NaN,'Color','b','linestyle','-','linewidth',2);

h(4) = plot(ax,NaN,NaN,'Color','r','linestyle','-','linewidth',2);


ylim(ax,[-0.05 1.05])
xlim(ax,[0 plotLen])

xlabel(ax,'Trial number')
ylabel(ax,'Fraction rightward choice')

legend(h,'Mouse','Model','Left DA','Right DA',...
   'location','northoutside','orientation','horizontal')


end