function Plot3BestWorst(Data,data_modelMin,data_modelMax)
% Find the 3 best and 3 worst Fvals and plot against

data = Data.data;


% changing y-axis values
data(:,3) = (1 + data(:,3)) ./ 2;
data_modelMin(:,:,16) = (1 + data_modelMin(:,:,16)) ./ 2;
data_modelMax(:,:,16) = (1 + data_modelMax(:,:,16)) ./ 2;


data_model_min1(:,:) = data_modelMin(1,:,:);
data_model_min2(:,:) = data_modelMin(2,:,:);
data_model_min3(:,:) = data_modelMin(3,:,:);
data_model_max1(:,:) = data_modelMax(1,:,:);
data_model_max2(:,:) = data_modelMax(2,:,:);
data_model_max3(:,:) = data_modelMax(3,:,:);


color=[1 0 0; 0 0 1; 0 1 0; 0 0 0];   % red, blue, green, black


% Initialise
per_mice = nan(1,length(unique(data(:,2))));
per_model_min1 = nan(1,length(unique(data_model_min1(:,2))));
per_model_min2 = nan(1,length(unique(data_model_min2(:,2))));
per_model_min3 = nan(1,length(unique(data_model_min3(:,2))));
per_model_max1 = nan(1,length(unique(data_model_max1(:,2))));
per_model_max2 = nan(1,length(unique(data_model_max2(:,2))));
per_model_max3 = nan(1,length(unique(data_model_max3(:,2))));

best = figure;
bestMice = subplot(1,4,1,'Parent',best); hold(bestMice,'on')
bestModel1 = subplot(1,4,2,'Parent',best); hold(bestModel1,'on')
bestModel2 = subplot(1,4,3,'Parent',best); hold(bestModel2,'on')
bestModel3 = subplot(1,4,4,'Parent',best); hold(bestModel3,'on')

worst = figure;
worstMice = subplot(1,4,1,'Parent',worst); hold(worstMice,'on')
worstModel1 = subplot(1,4,2,'Parent',worst); hold(worstModel1,'on')
worstModel2 = subplot(1,4,3,'Parent',worst); hold(worstModel2,'on')
worstModel3 = subplot(1,4,4,'Parent',worst); hold(worstModel3,'on')

for iblock=unique(data(:,8))'
   c=1;
   for istim = unique(data(:,2))'
      
      per_mice(c)=    nanmean(data(data(:,2)==istim & data(:,8)==iblock,3));
      
      per_model_min1(c)=   nanmean(data_model_min1(data_model_min1(:,2)==istim & data_model_min1(:,8)==iblock,16));
      per_model_min2(c)=   nanmean(data_model_min2(data_model_min2(:,2)==istim & data_model_min2(:,8)==iblock,16));
      per_model_min3(c)=   nanmean(data_model_min3(data_model_min3(:,2)==istim & data_model_min3(:,8)==iblock,16));
      
      per_model_max1(c)=   nanmean(data_model_max1(data_model_max1(:,2)==istim & data_model_max1(:,8)==iblock,16));
      per_model_max2(c)=   nanmean(data_model_max2(data_model_max2(:,2)==istim & data_model_max2(:,8)==iblock,16));
      per_model_max3(c)=   nanmean(data_model_max3(data_model_max3(:,2)==istim & data_model_max3(:,8)==iblock,16));
      
      c=c+1;
   end
   
   
   plot(bestMice,unique(data(:,2))',per_mice,'Color',color(iblock,:),...
      'marker','o','markersize',10,'markerfacecolor',color(iblock,:),...
      'linestyle','-','linewidth',1.2)
   
   plot(bestModel1,unique(data_model_min1(:,2))',per_model_min1,...
      'Color',color(iblock,:),'linestyle','--','linewidth',1.2,...
      'marker','o','markersize',10)
   
   plot(bestModel2,unique(data_model_min2(:,2))',per_model_min2,...
      'Color',color(iblock,:),'linestyle','--','linewidth',1.2,...
      'marker','o','markersize',10)
   
   plot(bestModel3,unique(data_model_min3(:,2))',per_model_min3,...
      'Color',color(iblock,:),'linestyle','--','linewidth',1.2,...
      'marker','o','markersize',10)
   
   
   plot(worstMice,unique(data(:,2))',per_mice,'Color',color(iblock,:),...
      'marker','o','markersize',10,'markerfacecolor',color(iblock,:),...
      'linestyle','-','linewidth',1.2)
   
   plot(worstModel1,unique(data_model_max1(:,2))',per_model_max1,...
      'Color',color(iblock,:),'linestyle','--','linewidth',1.2,...
      'marker','o','markersize',10)
   
   plot(worstModel2,unique(data_model_max2(:,2))',per_model_max2,...
      'Color',color(iblock,:),'linestyle','--','linewidth',1.2,...
      'marker','o','markersize',10)
   
   plot(worstModel3,unique(data_model_max3(:,2))',per_model_max3,...
      'Color',color(iblock,:),'linestyle','--','linewidth',1.2,...
      'marker','o','markersize',10)


end


   % plot titles
   title(bestMice,'Mice Raw Data')
   title(bestModel1,sprintf('Model with min Fval\n(Best)'))
   title(bestModel2,sprintf('Model with min2 Fval\n(2nd Best)'))
   title(bestModel3,sprintf('Model with min3 Fval\n(3rd Best)'))
   title(worstMice,'Mice Raw Data')
   title(worstModel1,sprintf('Model with max Fval\n(Worst)'))
   title(worstModel2,sprintf('Model with max2 Fval\n(2nd Worst)'))
   title(worstModel3,sprintf('Model with max3 Fval\n(3rd Worst)'))
   
   % plot y labels
   ylabel(bestMice,'Fraction rightward choice')
   ylabel(worstMice,'Fraction rightward choice')
   
   % plot x labels
   xlabel(bestMice,'Contrast')
   xlabel(bestModel1,'Contrast')
   xlabel(bestModel2,'Contrast')
   xlabel(bestModel3,'Contrast')
   xlabel(worstMice,'Contrast')
   xlabel(worstModel1,'Contrast')
   xlabel(worstModel2,'Contrast')
   xlabel(worstModel3,'Contrast')
   
   % plot y limits
   ylim(bestMice,[0 1])
   ylim(bestModel1,[0 1])
   ylim(bestModel2,[0 1])
   ylim(bestModel3,[0 1])
   ylim(worstMice,[0 1])
   ylim(worstModel1,[0 1])
   ylim(worstModel2,[0 1])
   ylim(worstModel3,[0 1])
   
   
end
