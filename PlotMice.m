function PlotMice(Data)

data = Data.data;

color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black

data(:,3) = (1 + data(:,3)) ./ 2;

per_mice = nan(1,length(unique(data(:,2))));

% PLOT
figure; hold on

for iblock = unique(data(:,8))'
   c=1;
   for istim = unique(data(:,2))'
      
      per_mice(c)= nanmean(data(data(:,2)==istim & data(:,8)==iblock,3)) ;
      
      c=c+1;
      
   end
   
   
   plot(unique(data(:,2))',per_mice,'Color',color(iblock,:),...
      'marker','o','markersize',10,'markerfacecolor','none',...
      'linestyle','-','color',color(iblock,:),'linewidth',1.2)
   ylabel('Fraction rightward choice')
   xlabel('Contrast')
   ylim([0 1])
   title ('Raw Mice Data')
   
end

legend('Block 1','Block 2','Block 3','Block 4','Location','southeast')

end
