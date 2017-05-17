function PlotModel(data_model)


color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black

% change the y-axis from (-1,1) to (0,1)
data_model(:,16) = (1 + data_model(:,16)) ./ 2;

per_model = nan(1,length(unique(data_model(:,2))));

% PLOT
figure; hold on

for iblock = unique(data_model(:,8))'
   c=1;
   for istim = unique(data_model(:,2))'
      
      per_model(c)= nanmean(data_model(data_model(:,2)==istim &...
         data_model(:,8)==iblock,16));
      
      c=c+1;
   end

   % plot alphaDiff = 0.0 in black
   plot(unique(data_model(:,2))',per_model,'Color',color(iblock,:),...
      'marker','o','markersize',10,'linestyle','-','linewidth',1.2)
   ylabel('Fraction rightward choice')
   xlabel('Contrast')
   ylim([0 1])
   title ('POMDP Model')
   
end

legend('Block 1','Block 2','Block 3','Block 4','Location','southeast')

end
