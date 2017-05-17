function PlotTrials(Data,data_model)
% nice plot

data_mice = Data.data;

% changing y-axis values
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
data_model(:,22) = (1 + data_model(:,22)) ./ 2;


plotLen = length(data_mice);
if length(data_mice)<plotLen
   
   plotLen = length(data_mice);
   
end
plotLenVec = cat(1,ones(plotLen,1),zeros(length(data_mice)-plotLen,1));

mouse = data_mice(plotLenVec==1,3) ;
model = data_model(plotLenVec==1,22);

smoothmouse = mouse*0.96 + 0.02;
smoothmodel = model*0.96 + 0.02;

smoothmouse = smooth(smoothmouse,7);
smoothmodel = smooth(smoothmodel,7);


blockID = data_mice(plotLenVec==1,8) - 1;
blockID(blockID==1) = 1.03;
blockID(blockID==0) = -0.03;

figure;
subplot(2,1,1);  hold on
plot(mouse,'marker','o','MarkerEdgeColor',[0.2 0 0.6],'markerfacecolor',...
   [0.2 0 0.6],'markersize',1.5,'linestyle','none')

plot(blockID,'marker','o','MarkerEdgeColor',[1 .5 0],'markerfacecolor',...
   [1 .5 0],'markersize',2,'linestyle','none')

plot(smoothmouse,'color',[0.2 0 0.6])

ylim([-0.05 1.05])
xlim([0 plotLen])

xlabel('Trial number')
ylabel('Fraction rightward choice')
title('Trial by trial mouse actions')


subplot(2,1,2); hold on
plot(mouse,'marker','o','MarkerEdgeColor',[0.2 0 0.6],'markerfacecolor',...
   [0.2 0 0.6],'markersize',1.5,'linestyle','none')

plot(blockID,'marker','o','MarkerEdgeColor',[1 .5 0],'markerfacecolor',...
   [1 .5 0],'markersize',2,'linestyle','none')

%plot(smoothmodel,'color',[0.7 0.7 0.7])
plot(smoothmodel,'color',[0 0 0])

ylim([-0.05 1.05])
xlim([0 plotLen])

xlabel('Trial number')
ylabel('Fraction rightward choice')
title('Trial by trial model actions')

set(gca,'visible','on');


end