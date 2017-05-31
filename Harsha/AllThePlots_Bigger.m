 function AllThePlots_Bigger(Data, data_modelAll, actionAll, correctAll, xanswerMin, includeContrast, NLLData, varargin)

 % 1 - (Left) Psychometric curves of blocks 1 and 2. (Right) Details of parameters fit
 % 2 -  Model and experimental P(R) sliding over all sessions
 % 3 -  NLL : (Left) Alpha v/s Noise  (Middle) DA v/s Noise  (Right) Bias v/s Noise
 % 
 % If blocks 3/4 present:
 % 4 - (Left/centre) Conditional psychometric curves of Block 3/4 (Right) Alpha with DA conditional NLL
 
 
 % 4/5 - (Left) Percent match for top 5 models  (Right) Percent shift in P(R) of expt and top 5 models
 % Combining all the plots we have

data_mice = Data.data;
animalID = Data.ID;

% changing y-axis values
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
data_modelAll(:,16,:) = (1 + data_model(:,16,:)) ./ 2;
data_modelAll(:,22,:) = (1 + data_model(:,22,:)) ./ 2;
actionAll = ( 1 + actionAll) ./2;

% Best model
data_model = squeeze(data_modelAll(:,:,1));      
action = actionAll(:,:,1);
correct = correctAll(:,:,1);

nModels = size(data_modelAll, 3);
nTrials = size(data_modelAll,1);
nIters = size(action,2);

% Plotting and data parameters
color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black
markerShape = {'o', 's', '^', 'v'};
BlockLabels = {'DA-L', 'DA-R','DA-Both', 'No DA'}; 

if any(data_mice(:,8)>2)
    % Block 3/4 present
    nrows = 5;
else
    nrows = 4;
end

contrasts = unique(data_mice(:,2))';
if nargin<8    
    AllBlocks = unique(data_mice(:,8))';
else
    AllBlocks = varargin{1};
end

includeContrast = reshape( includeContrast, length(contrasts), length(AllBlocks) );

FvalStore   = NLLData.FvalStore;
alpha       = NLLData.alpha;
DAval       = NLLData.DAval;
noiseSTD    = NLLData.noiseSTD;
Qbias       = NLLData.Qbias;

alpha2      = NLLData.alpha2;
FvalStore2  = NLLData.alpha2;
params2     = NLLData.top5params;

full = figure;

%-------------------------------------------------------------------------%
%% PLOT 1 - (Left) Psychometric curves of blocks 1 and 2. (Right) Details of parameters fit
%-------------------------------------------------------------------------%

%-------------- Left -----------------------------------------------------%
subplot(nrows, 2, 1 ); hold on;
blocks = [ 1 2];
b=0;
h = zeros(length(blocks)+2,1);
per_mice = nan(1, length(contrasts) );
per_model = per_mice;
for blockID = blocks
   b=b+1;
   c=1;
   for istim = contrasts
      
      id=(data_mice(:,2)==istim & data_mice(:,8)==blockID);
      per_mice(c)= nanmean(data_mice(id,3)) ;
      
      per_model(c)=0;
      for iter = 1:nIters
          id =(data_model(:,2)==istim & data_model(:,8)==blockID);
          per_model(c)= per_model(c) + nanmean(action( id,iter)) ;
      end
      per_model(c)=per_model(c)/nIters;
      c=c+1;
   end
   
   plot_contrast = includeContrast(:, blockID )==1;
  
   
   % PLOT
   % plot the mice raw data as points
  
   plot(contrasts(plot_contrast),per_mice(plot_contrast),...
      'color',color(blockID,:),...
      'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color(blockID,:),...
      'markerfacecolor',color(blockID,:),'linestyle','-',...
      'linewidth',1.2);
   plot(contrasts(plot_contrast),per_model(plot_contrast),'color',color(blockID,:),...
      'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color(blockID,:),...
      'linestyle','--','linewidth',1.2);

   h( b ) = plot(ax,NaN,NaN,'color',color(blockID,:),'marker',markerShape{blockID},'markersize',10,...
  'markerfacecolor',color(blockID,:),'markerfacecolor',color(blockID,:),'linestyle','none');
   

end
   
% dummy plots for the legend
b=length(blocks);
h(b+1) = plot(ax,NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
h(b+2) = plot(ax,NaN,NaN,'color','k','linestyle','--','linewidth',1.2);

legend(h, BlockLabels{blocks},'Mouse','Model','Location','southeast');
ylabel(ax,'Fraction rightward choice')
xlabel(ax,'Contrast')
title('Psychometric curves from Data and Best Model')

%---------------------- Right --------------------------------------------%
% put all the explanation in this subplot
subplot(nrows,2,2); hold on;
set(gca,'visible','off');

descr1 = {['Animal ID:  ',animalID]};
text(0,0.9,descr1)

descr2 = {'Model parameter values:';
   ['alpha = ',num2str(xanswerMin(1,1))];
   ['DAval = ',num2str(xanswerMin(1,2))];
   ['noiseSTD = ',num2str(xanswerMin(1,3))];
   ['Bias in QValue = ',num2str(xanswerMin(1,4))];};
text(.55,0.9,descr2)

%-------------------------------------------------------------------------%
%% PLOT 2 -  Model and experimental P(R) sliding over all sessions
%-------------------------------------------------------------------------%


subplot(nrows,1,2,'Parent',full);  hold on

plotLen = length(data_mice);
mouse = data_mice(:,3) ;
model = data_model(:,22);

% Y Lim of choices put between 0.02 and 0.98
smoothmouse = mouse*0.96 + 0.02;
smoothmodel = model*0.96 + 0.02;

% Smoothed over 15 trials
smoothmouse = smooth(smoothmouse,15);
smoothmodel = smooth(smoothmodel,15);


blockID = data_mice(:,8);
mouse_color = 'k';
model_color = [0.5 0 0.5];

%X and O at 0/1 for mouse and model
plot(mouse,'marker','x','MarkerEdgeColor',mouse_color,'markerfacecolor',...
   mouse_color,'markersize',2,'linestyle','none')

plot(model,'marker','o','MarkerEdgeColor',model_color,...
    'markersize',2,'linestyle','none')

% Block Identities: blue for DA-L, red for DA-R, green for DA-both
plot( find(blockID==1),-0.03, 'marker', 's', 'markeredgecolor', 'b', 'markerfacecolor', 'b','linestyle', none );
plot( find(blockID==2), 1.03, 'marker', 's', 'markeredgecolor', 'r', 'markerfacecolor', 'r','linestyle', none );
plot( find(blockID==3), 1.03, 'marker', 's', 'markeredgecolor', 'g', 'markerfacecolor', 'g','linestyle', none );
plot( find(blockID==3),-0.03, 'marker', 's', 'markeredgecolor', 'g', 'markerfacecolor', 'g','linestyle', none );

% Plot smoothed P(R)
plot(smoothmouse,'color',mouse_color)
plot(smoothmodel,'color',model_color)


ylim([-0.05 1.05])
xlim([0 plotLen])

xlabel('Trial number')
ylabel('Fraction rightward choice')

set(gca,'visible','on');

%-------------------------------------------------------------------------%
%% PLOT 3 - NLL : (Left) Alpha v/s Noise  (Middle) DA v/s Noise  (Right) Bias v/s Noise
%-------------------------------------------------------------------------%

subplot(nrows,1,3)
title('NLL of Psychometric curves')
[~,locMin] = min(FvalStore(:));
[a,b,c,d] = ind2sub(size(FvalStore),locMin);

%------------- (Left) Alpha v/s Noise ------------------%
clear colourVal x y
colourVal = nan( length(alpha),length(noiseSTD));
for j=1:length(noiseSTD)
   for i=1:length(alpha)
      
      % Set labels for the parameters we are varying
      xData = 'noiseSTD';
      yData = 'alpha';
      
      % Prepare variables for colour map
      x = noiseSTD; %[min(noiseSTD),max(noiseSTD)];
      y = alpha;    %[min(alpha),max(alpha)];
      colourVal(i,j) = FvalStore(i,b,j,d);
      
   end
end

subplot(nrows,3,7);
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)

%------------- (Middle) DAVal v/s Noise ------------------%
clear colourVal x y
colourVal = nan(length(DAval),length(noiseSTD));
for j=1:length(noiseSTD)
   for i=1:length(DAval)
      
      % Set labels for the parameters we are varying
      xData = 'noiseSTD';
      yData = 'DAval';
      
      % Prepare variables for colour map
      x = noiseSTD; %[min(noiseSTD),max(noiseSTD)];
      y = DAval;    %[min(DAval),max(DAval)];
      colourVal(i,j) = FvalStore(a,i,j,d);
      
   end
end

% plot the colour map
subplot(nrows,3,8)
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)

%------------- (Right) Qbias v/s Noise ------------------%
clear colourVal x y
colourVal = nan(length(Qbias),length(noiseSTD));
for j=1:length(noiseSTD)
   for i=1:length(Qbias)
      
      % Set labels for the parameters we are varying
      xData = 'noiseSTD';
      yData = 'Qbias';
      
      % Prepare variables for colour map
      x = noiseSTD; %[min(noiseSTD),max(noiseSTD)];
      y = Qbias;    %[min(Qbias),max(Qbias)];
      colourVal(i,j) = FvalStore(a,b,j,i);
      
   end
end

% plot the colour map
subplot(nrows,3,9)
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)

%-------------------------------------------------------------------------%
%% (Conditional) PLOT 4 - Conditional choices (For Blocks 3/4 if PRESENT)
%-------------------------------------------------------------------------%

%-- (Left) Conditional Psychometric curves after Left and Right choices
subplot( nrows, 2, 7); hold on;
if any(data_mice(:,8)>2)
   % Block 3/4 present
   blocks = unique(data_mice(data_mice(:,8)>2,8))'; 
end

per_mice = nan(1, length(contrasts) );
per_model = per_mice;

h = zeros(2*length(blocks)+2,1);
b=0;
color_b={'b', 'cyan'};
color_r = {'r', 'm'};
hLabel = cell(2*length(blocks),1);
for blockID = blocks
   b=b+1;
   c=1;
   for istim = contrasts
      
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
      
      c=c+1;
   end

   plot_contrast = includeContrast(:, blockID )==1;
   
   % PLOT
   plot(contrasts(plot_contrast),per_mice(1,plot_contrast),...
      'color',color_b{b},...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor',color_b{b},...
      'markerfacecolor',color_b{b},'linestyle','-',...
      'linewidth',1.2);
   plot(ax,contrasts(plot_contrast),per_mice(2,plot_contrast),...
      'color',color_r{b},...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor',color_r{b},...
      'markerfacecolor',color_r{b},'linestyle','-',...
      'linewidth',1.2);

   
   % plot the model
   plot(ax,contrasts(plot_contrast),per_model(1,plot_contrast),'color',color_b{b},...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor',color_b{b},...
      'linestyle','--','linewidth',1.2);
   plot(ax,contrasts(plot_contrast),per_model(2,plot_contrast),'color',color_r{b},...
      'marker',markerShapes{blockID},'markersize',10,'markeredgecolor',color_r{b},...
      'linestyle','--','linewidth',1.2);

   ylabel(ax,'Fraction rightward choice')
   xlabel(ax,'Contrast')
   h( 2*b-1 ) = plot(NaN,NaN,'color',color_b{b},'marker',markerShape{blockID},'markersize',10,...
   'markerfacecolor',color_b{b},'markerfacecolor',color_b{b},'linestyle','-');
   h( 2*b ) = plot(NaN,NaN,'color',color_r{b},'marker',markerShape{blockID},'markersize',10,...
   'markerfacecolor',color_r{b},'markerfacecolor',color_r{b},'linestyle','-');
  
   hLabel{2*b-1} = strcat( 'After Left - ', BlockLabels{blockID});
   hLabel{2*b}   = strcat( 'After Right - ', BlockLabels{blockID});
end

   b = 2*length(blocks);
   
 
   h(b+1) = plot(ax,NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
   h(b+2) = plot(ax,NaN,NaN,'color','k','linestyle','--','linewidth',1.2);

   legend(h, hLabel{:}, 'Mouse','Model','Location','southeast');
% end


%-- (Right) NLL of conditional PC as a function of alpha
subplot( nrows, 2, 8); hold on;
plotDA = params2(:,2);
imagesc( plotDA, alpha2, FvalStore2 );
title(' NLL of conditional PC for top 5 models (Columns)');

end

