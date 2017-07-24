 function AllThePlots_Bigger_WithLapse(Data, data_modelAll, actionAll, correctAll, xanswerMin, includeContrast, NLLData, varargin)

 % 1 - (Left) Psychometric curves of blocks 1 and 2. (Middle) Percent match for top 5 models (Right) Details of parameters fit
 % 2 -  Model and experimental P(R) sliding over all sessions
 % 3 -  NLL : (Left) Alpha v/s Noise  (Middle) DA v/s Noise  (Right) Bias v/s Noise
 % 
 % If blocks 3/4 present:
 % 4 - (Left/centre) Conditional psychometric curves of Block 3/4 (Right) Alpha with DA conditional NLL
 % 
 % 4/5 - (Left) Shifts after easy/difficult stimuli - mouse and model. (Right) Sifts after DA/water reward and error - mouse and model
 % Combining all the plots we have

data_mice = Data.data;
animalID  = Data.ID;

% changing ACTION values to 0 and 1
if any(data_mice(:,3)==-1)
data_mice(:,3)        = (1 + data_mice(:,3)) ./ 2;
end
if any(data_modelAll(:,22,1)==-1)
data_modelAll(:,16,:) = (1 + data_modelAll(:,16,:)) ./ 2;
data_modelAll(:,22,:) = (1 + data_modelAll(:,22,:)) ./ 2;
end
if any(actionAll(:,1,1)==-1)
actionAll = ( 1 + actionAll) ./2;
end

% Best model
data_model = squeeze(data_modelAll(:,:,1));      
action     = actionAll(:,:,1);
correct    = correctAll(:,:,1);

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
lambdaL     = NLLData.lambdaL;
lambdaR     = NLLData.lambdaR;


alpha2      = NLLData.alpha2;
FvalStore2  = NLLData.FvalStore2;
params2     = NLLData.top5params;

full = figure;

%-------------------------------------------------------------------------%
%% PLOT 1 - MODEL DETAILS  
%-------------------------------------------------------------------------%

%% -------- (Left) Psychometric curves of blocks 1 and 2 ----------------%
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

   h( b ) = plot(NaN,NaN,'color',color(blockID,:),'marker',markerShape{blockID},'markersize',10,...
  'markerfacecolor',color(blockID,:),'markerfacecolor',color(blockID,:),'linestyle','none');
   

end
   
% dummy plots for the legend
b=length(blocks);
h(b+1) = plot(NaN,NaN,'color','k','linestyle','-','linewidth',1.2);
h(b+2) = plot(NaN,NaN,'color','k','linestyle','--','linewidth',1.2);

legend(h, BlockLabels{blocks},'Mouse','Model','Location','southeast');
ylabel('Fraction rightward choice')
ylim([0 1])
xlabel('Contrast')
title('Psychometric curves from Data and Best Model')

%% ----- (Middle) Percent match for top 5 models --------------%
PMatch = nan(nIters, nModels);
PMatch2 = nan(nIters, nModels);     % match on non-zero contrasts
PMatch_mode = nan(nModels,1);
id = (data_mice(:,2) ~= 0  );
nCTrials = sum(id);
for model = 1:nModels
    for jj=1:nIters
      PMatch(jj, model)  = sum( data_mice(:,3) == actionAll(:, jj,model) )/nTrials; 
      PMatch2(jj, model) = sum( data_mice(id,3) == actionAll(id, jj,model) )/nCTrials; 

    end
    PMatch_mode(model) = sum(data_mice(:,3) == data_modelAll(:,16,model))/nTrials;
end

subplot(nrows,4,3); hold on

h = zeros(2,1);
%%% only plotting mean percent match for ALL trials
h(1) = plot(1:nModels, mean(PMatch,1), 'color', 'r', 'marker', 'o', 'markerfacecolor','r', 'markersize', 10,'linestyle','none');
% h(2) = plot(1:nModels, mean(PMatch,1), 'color', 'cyan', 'marker', 'o', 'markerfacecolor','cyan', 'markersize', 10,'linestyle','none');
% h(2) = plot(1:nModels, PMatch_mode, 'color', 'cyan', 'marker', 'o', 'markerfacecolor','cyan', 'markersize', 10);
for jj=1:nIters
    plot(1:nModels, PMatch(jj,:), 'marker', '*', 'markerfacecolor', 'k', 'linestyle','none', 'markersize',4);
%     plot(1:nModels, PMatch2(jj,:),'marker', 's', 'markerfacecolor', 'k', 'linestyle','none', 'markersize',4);
end
legend(h(1), 'Mean Percent Match of responses');
% legend(h, 'Mean Percent Match of responses', 'Percent response match for non-zero contrast trials')%, 'Percent match of mode response')
set(gca, 'xtick', 1:nModels);
xlabel( 'Model number')

%% ------------ (Right) Details of parameters fit ------------------------%
% put all the explanation in this subplot
subplot(nrows,4,4); hold on;
set(gca,'visible','off');

descr1 = {['Animal ID:  ',animalID]};
text(0,0.9,descr1)

descr2 = {'Model parameter values:';
   ['alpha = ',num2str(xanswerMin(1,1))];
   ['DAval = ',num2str(xanswerMin(1,2))];
   ['noiseSTD = ',num2str(xanswerMin(1,3))];
   ['Bias in QValue = ',num2str(xanswerMin(1,4))];
   ['LambdaL = ',num2str(xanswerMin(1,5))];
   ['LambdaR = ',num2str(xanswerMin(1,6))];};
text(.55,0.9,descr2)



%-------------------------------------------------------------------------%
%% PLOT 2 -  Model and experimental P(R) sliding over all sessions
%-------------------------------------------------------------------------%

subplot(nrows,1,2,'Parent',full);  hold on

plotLen = nTrials;
mouse   = data_mice(:,3) ;
model   = data_model(:,22);

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
plot( find(blockID==1),-0.03, 'marker', 's', 'markeredgecolor', 'b', 'markerfacecolor', 'b','linestyle', 'none' ,'markersize',2);
plot( find(blockID==2), 1.03, 'marker', 's', 'markeredgecolor', 'r', 'markerfacecolor', 'r','linestyle', 'none' ,'markersize',2);
if any(blockID==3)
    plot( find(blockID==3), 1.03, 'marker', 's', 'markeredgecolor', 'g', 'markerfacecolor', 'g','linestyle', 'none' ,'markersize',2);
plot( find(blockID==3),-0.03, 'marker', 's', 'markeredgecolor', 'g', 'markerfacecolor', 'g','linestyle', 'none' ,'markersize',2);
end

% Plot smoothed P(R)
plot(smoothmouse,'color',mouse_color)
plot(smoothmodel,'color',model_color)


ylim([-0.05 1.05])
xlim([0 plotLen])

xlabel('Trial number')
ylabel('Fraction\newlinerightward choice')

set(gca,'visible','on');

%-------------------------------------------------------------------------%
%% PLOT 3 - NLL of conditional Psychometric curves
%-------------------------------------------------------------------------%

subplot(nrows,1,3)
title('NLL of Psychometric curves')
[~,locMin] = min(FvalStore(:));
[a,b,~,d,e,f] = ind2sub(size(FvalStore),locMin);

%% ------------- (Left) Alpha v/s Noise ------------------%
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
      colourVal(i,j) = FvalStore(i,b,j,d,e,f);
      
   end
end

subplot(nrows,3,7);
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)

%% ------------- (Middle) DAVal v/s Noise ------------------%
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
      colourVal(i,j) = FvalStore(a,i,j,d,e,f);
      
   end
end

% plot the colour map
subplot(nrows,3,8)
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)

%% ------------- (Right) Qbias v/s Noise ------------------%
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
      colourVal(i,j) = FvalStore(a,b,j,i,e,f);
      
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
if any(data_mice(:,8)>2)
%-- (Left) Conditional Psychometric curves after Left and Right choices
   subplot( nrows, 2, 7); hold on;

   % Block 3/4 present
   blocks = unique(data_mice(data_mice(:,8)>2,8))'; 
   mB = 5;
   if ~isempty(blocks)
       mB = max(blocks);
   end

per_mice = nan(1, length(contrasts) );
per_model = per_mice;

% h = zeros(2*length(blocks)+2,1);
h = zeros(2+2,1);
b=0;
color_b={'cyan','cyan'};
color_r = {'m','m'};
% hLabel = cell(2*length(blocks),1);
hLabel = cell(2,1);
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
      
      for iter = 1:nIters
          id =[false;data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==0 & correct(1:end-1,iter)==1];
          per_model(1,c)= per_model(1,c) + nanmean(action( id,iter)) ;
          id = [false; data_model(2:end,2)==istim & data_model(2:end,8)==blockID & ...
             action(1:end-1,iter)==1 & correct(1:end-1,iter)==1];
          per_model(2,c)= per_model(2,c) + nanmean(action( id,iter)) ;
      end
      per_model(1,c)=per_model(1,c)/nIters;
      per_model(2,c)=per_model(2,c)/nIters;
      
      c=c+1;
   end

   plot_contrast = includeContrast(:, 2+b )==1;
   
   % PLOT only block 4, or block 3 if no block 4 present
   if blockID == mB
       
       % plot mouse
       plot(contrasts(plot_contrast),per_mice(1,plot_contrast),...
          'color',color_b{b},...
          'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color_b{b},...
          'markerfacecolor',color_b{b},'linestyle','-',...
          'linewidth',1.2);
       plot(contrasts(plot_contrast),per_mice(2,plot_contrast),...
          'color',color_r{b},...
          'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color_r{b},...
          'markerfacecolor',color_r{b},'linestyle','-',...
          'linewidth',1.2);


       % plot the model
       plot(contrasts(plot_contrast),per_model(1,plot_contrast),'color',color_b{b},...
          'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color_b{b},...
          'linestyle','--','linewidth',1.2);
       plot(contrasts(plot_contrast),per_model(2,plot_contrast),'color',color_r{b},...
          'marker',markerShape{blockID},'markersize',10,'markeredgecolor',color_r{b},...
          'linestyle','--','linewidth',1.2);

       ylabel('Fraction\newlinerightward choice')
       xlabel('Contrast')
       h( 1 ) = plot(NaN,NaN,'color',color_b{b},'marker',markerShape{blockID},'markersize',10,...
       'markerfacecolor',color_b{b},'markerfacecolor',color_b{b},'linestyle','-');
       h( 2 ) = plot(NaN,NaN,'color',color_r{b},'marker',markerShape{blockID},'markersize',10,...
       'markerfacecolor',color_r{b},'markerfacecolor',color_r{b},'linestyle','-');

       hLabel{1} = strcat( 'After Left - ', BlockLabels{blockID});
       hLabel{2}   = strcat( 'After Right - ', BlockLabels{blockID});
   end
end

   %b = 2*length(blocks);
   
 
   h(3) = plot(NaN,NaN,'color','k','linestyle','-','linewidth',1.2);    
   h(4) = plot(NaN,NaN,'color','k','linestyle','--','linewidth',1.2);

   legend(h, hLabel{:}, 'Mouse','Model','Location','southeast');
% end


%% -- (Right) NLL of conditional PC as a function of alpha
subplot( nrows, 2, 8); hold on;
% plotDA = params2(:,2);
% imagesc( plotDA, alpha2, FvalStore2 );


x = params2(:,2)'; 
y = alpha2; 
xData = 'DA value of model';
yData = 'Alpha';

% plot the colour map
imagesc(x,y,FvalStore2);
colorbar
xlabel(xData)
ylabel(yData)
title(' NLL of conditional PC for top 5 models (Columns)');
end


%-------------------------------------------------------------------------%
%% PLOT 5 - Performance of top models  
%-------------------------------------------------------------------------%
%data_model(:,16) --> mode response/ data_model(:,22) --> Percent match?



%% ------------ (Left) Single Trial Shift (STS) in P(R) after dopamine reward : Conditioned on stimulus difficulty ----%

model = 1; % only compare best model
% shift
[ stimuli_M, mouse_shift_allstim, mouse_ED_PC, mouse_ED_IC ] = SingleTrialShift_EasyDiffStim_with_Dop( data_mice );
[ stimuli_m, model_shift_allstim, model_ED_PC, model_ED_IC ] = SingleTrialShift_EasyDiffStim_with_Dop( squeeze(data_modelAll(:,:,model)), squeeze(actionAll(:,:,model)), squeeze(correctAll(:,:,model))  );
PlotShift = max(stimuli_M);


subplot(nrows,2,2*nrows-1); hold on
h=zeros(8+1,1);
MLabel = cell(8+1,1);
% PLOT MOUSE
% delta P(R) will be negative
h(1) = plot( stimuli_M, mouse_shift_allstim(:,1,1),  'color', 'k', 'marker', 'o', 'markerfacecolor', 'k', 'markersize', 5, 'linestyle', ':',  'linewidth',2) ;
h(2) = plot( stimuli_M, mouse_shift_allstim(:,1,2),  'color', 'b', 'marker', 'o', 'markerfacecolor', 'b', 'markersize', 5, 'linestyle', '-',  'linewidth',2) ;
% delta P(R) will be positive
h(3) = plot( stimuli_M, mouse_shift_allstim(:,2,1),  'color', 'k', 'marker', 'o', 'markerfacecolor', 'k', 'markersize', 5, 'linestyle', '--', 'linewidth',2) ;
h(4) = plot( stimuli_M, mouse_shift_allstim(:,2,2),  'color', 'r', 'marker', 'o', 'markerfacecolor', 'r', 'markersize', 5, 'linestyle', '-', 'linewidth',2) ;
MLabel{1} = 'Mouse - Easy Stim in DA-L';
MLabel{2} = 'Diff stim in DA-L';
MLabel{3} = 'Easy stim in DA-R';
MLabel{4} = 'Diff stim in DA-R';

h(5) = plot( NaN, NaN, 'color', 'w'); %Blank  line in legend
MLabel{5} = '';

% PLOT BEST MODEL
model = 1; % only plot model 1 
h(6) = plot( 1+PlotShift + stimuli_m, nanmean(model_shift_allstim(:,1,1,:),4),  'color', 0.5*[1 1 1], 'marker', 's',  'markerfacecolor', 0.5*[1 1 1],'markersize', 5, 'linestyle', ':',  'linewidth',2) ;
h(7) = plot( 1+PlotShift + stimuli_m, nanmean(model_shift_allstim(:,1,2,:),4),  'color', 'cyan',      'marker', 's',  'markerfacecolor', 'cyan',     'markersize', 5, 'linestyle', '-',  'linewidth',2) ;
                
h(8) = plot( 1+PlotShift + stimuli_m, nanmean(model_shift_allstim(:,2,1,:),4),  'color', 0.5*[1 1 1], 'marker', 's',  'markerfacecolor', 0.5*[1 1 1],'markersize', 5, 'linestyle', '--', 'linewidth',2) ;
h(9) = plot( 1+PlotShift + stimuli_m, nanmean(model_shift_allstim(:,2,2,:),4),  'color', 'm',         'marker', 's',  'markerfacecolor', 'm',        'markersize', 5, 'linestyle', '-', 'linewidth',2) ;
                
MLabel{6} = sprintf('Model %d - Easy Stim in DA-L', model);
MLabel{7} = 'Diff stim in DA-L';
MLabel{8} = 'Easy stim in DA-R';
MLabel{9} = 'Diff stim in DA-R';

xlim([min(stimuli_m)-0.2 1+PlotShift+ max(stimuli_m)+PlotShift])
ylim( [-max(abs(mouse_shift_allstim(:)))-0.1 max(abs(mouse_shift_allstim(:)))+0.1])
ylabel('Change in P(R)')
xlabel('Stimulus')
legend(h, MLabel{:}, 'Location', 'NorthEast');
set(gca, 'xtick', ([stimuli_M(2:2:end), 1+PlotShift+stimuli_M(2:2:end)]), 'xticklabel',string([stimuli_M(2:2:end), stimuli_M(2:2:end)]) )

%% ------------ (Right) Single Trial Shift (STS) after Dopamine or Water Reward and Omission ----%

% on dopamine/water reward trials, or error trials of actions associated with
% dopamine/water

% new independent way of calculating shift - doesn't have stim levels yet
% mouse_shift_dop = CondPsychCurves_AfterDop( data_mice(:,2), data_mice(:,8), data_mice(:,3), data_mice(:,10) );
% mouse_shift_water = CondPsychCurves_AfterWater( data_mice(:,2), data_mice(:,8), data_mice(:,3), data_mice(:,10) );

%shift : [ stimulus level -2 to 2,     correct/error,      blockID,       numIters ]
mouse_shift_dop     = SingleTrialShift_DopTrials_OneEasy( data_mice );
mouse_shift_water   = SingleTrialShift_WaterTrials_OneEasy( data_mice );


% delta P(R) on block 2 + delta P(L) on block 1
mouse_dopreward_STS     = squeeze(mouse_shift_dop(:,1,2,:) - mouse_shift_dop(:,1,1,:));               % expected to be positive
mouse_doperror_STS      = squeeze(mouse_shift_dop(:,2,2,:) - mouse_shift_dop(:,2,1,:));               % expected to be negative

% delta P(L) on block 2 + delta P(R) on block 1
mouse_waterreward_STS   = squeeze( - mouse_shift_water(:,1,2,:) + mouse_shift_water(:,1,1,:));        % expected to be positive
mouse_watererror_STS    = squeeze( - mouse_shift_water(:,2,2,:) + mouse_shift_water(:,2,1,:));        % expected to be negative

%shift : [ stimulus level -2 to 2,     correct/error,      blockID,       numIters,         numModels ]
model_shift_dop = nan( [5, 2, 2, 99, 5] );   %Shifts for each top model and its iterations  
model_shift_water = model_shift_dop;

for model = 1:5
  model_shift_dop(:,:,:,:,model)   = SingleTrialShift_DopTrials_OneEasy( squeeze(data_modelAll(:,:,model)), ...
                                        squeeze(actionAll(:,:,model)), squeeze(correctAll(:,:,model))  );
  model_shift_water(:,:,:,:,model) = SingleTrialShift_WaterTrials_OneEasy( squeeze(data_modelAll(:,:,model)), ...
                                        squeeze(actionAll(:,:,model)), squeeze(correctAll(:,:,model))  );
                                    
  % delta PR on block 2 + delta PL on block 1
  model_dopreward_STS   = squeeze(   model_shift_dop(:,1,2,:,:)   - model_shift_dop(:,1,1,:,:));        % expected to be positive
  model_doperror_STS    = squeeze(   model_shift_dop(:,2,2,:,:)   - model_shift_dop(:,2,1,:,:));        % expected to be negative
  
  % delta PL on block 2 + delta PR on block 1
  model_waterreward_STS = squeeze( - model_shift_water(:,1,2,:,:) + model_shift_water(:,1,1,:,:));       % expected to be positive
  model_watererror_STS  = squeeze( - model_shift_water(:,2,2,:,:) + model_shift_water(:,2,1,:,:));       % expected to be negative

end

% shifts combined across blocks: 
%model - [stimulus, iter, modelNumber], Mice - [stimulus]


subplot(nrows,2,2*nrows); hold on
h=zeros(8+1,1);
MLabel = cell(8+1,1);

%%% PLOT MOUSE
h(1) = plot( -2:2, mouse_dopreward_STS,                            'color', [0 0.5 0],   'marker', 'o', 'markerfacecolor', [0 0.5 0], 'markersize', 10, 'linestyle', '-',  'linewidth',2) ;
h(2) = plot( -2:2, mouse_doperror_STS,                             'color', 'r',         'marker', 'o', 'markerfacecolor', 'r',       'markersize', 10, 'linestyle', '-',  'linewidth',2) ;
h(3) = plot( -2:2, mouse_waterreward_STS,                          'color', [0 0.5 0],   'marker', 'o',                               'markersize', 10, 'linestyle', '--', 'linewidth',2) ;
h(4) = plot( -2:2, mouse_watererror_STS,                           'color', 'r',         'marker', 'o',                               'markersize', 10, 'linestyle', '--', 'linewidth',2) ;
MLabel{1} = 'Mouse - DA reward';
MLabel{2} = 'DA Error';
MLabel{3} = 'Water reward';
MLabel{4} = 'Water Error';

h(5) = plot( NaN, NaN, 'color', 'w'); %Blank  line in legend
MLabel{5} = '';

%%% PLOT BEST MODEL
model = 1; % only plot model 1
h(6) = plot( 6+(-2:2), nanmean(model_dopreward_STS(:,:,model),2),  'color', 'g',         'marker', 's',  'markerfacecolor', 'g',   'markersize', 10, 'linestyle', '-',  'linewidth',1) ;
h(7) = plot( 6+(-2:2), nanmean(model_doperror_STS(:,:,model),2),   'color', 'm',         'marker', 's',  'markerfacecolor', 'm',   'markersize', 10, 'linestyle', '-',  'linewidth',1) ;                
h(8) = plot( 6+(-2:2), nanmean(model_waterreward_STS(:,:,model),2),'color', 'g',         'marker', 's',                            'markersize', 10, 'linestyle', '--', 'linewidth',1) ;
h(9) = plot( 6+(-2:2), nanmean(model_watererror_STS(:,:,model),2), 'color', 'm',         'marker', 's',                            'markersize', 10, 'linestyle', '--', 'linewidth',1) ;                
MLabel{6} = sprintf('Model %d - DA reward', model);
MLabel{7} = 'DA Error';
MLabel{8} = 'Water reward';
MLabel{9} = 'Water Error';
xlim([-3 3+7])
legend(h, MLabel{:}, 'Location', 'NorthEast');
set(gca, 'xtick', [-2:2 6+(-2:2)], 'xticklabel', {'Easy\newlineLeft', 'Diff\newlineLeft', 'Zero\newlineContrast', 'Diff\newlineRight', 'Easy\newlineRight', 'Easy\newlineLeft', 'Diff\newlineLeft', 'Zero\newlineContrast', 'Diff\newlineRight', 'Easy\newlineRight' } )
ylabel('Change in\newlineP[previous action]')
xlabel('Stimulus')

 end