
addpath(genpath('../../POMDP'));


close all


% set model parameter values
% alpha       = 0.0 : 0.05 : 1.0;      % learning rate
% DA_val      = 0.0 : 0.25 : 8.0;      % dopamine value
% noiseSTD    = 0.01 : 0.04 : 0.91;    % noise in belief

% TESTING2 set model parameter values
alpha       = 0.0 : 0.03 : 0.7;      % learning rate
DA_val      = 0.0 : 0.3 : 8.0;      % dopamine value
noiseSTD    = 0.02 : 0.02 : 0.42;    % noise in belief

% TESTING set model parameter values
% alpha       = 0.0 : 0.5 : 1.0;      % learning rate
% DA_val      = 0.0 : 3.0 : 6.0;      % dopamine value
% noiseSTD    = 0.03 : 0.2 : 0.43;    % noise in belief

% load in animal data
AnimalID='ALK028'; % 'SS031_DA' or 'SS040' or 'ALK05' or 'ALK011'
% or 'ALK028' or 'ALK017'

ExpID = '1';

data = LoadSavedDataForBehExp (AnimalID, ExpID);

% blockIDs 1 & 2 represent asymmetric dopamine reward
% remove blockIDs 3 & 4
data(data(:,8)==3,:)=[];
data(data(:,8)==4,:)=[];

% data(data(:,8)==1,:)=[];
% data(data(:,8)==2,:)=[];


% check number of trials for each contrast
contrast = unique(data(:,2))';
trialsPerContrast = nan(length(contrast),1);
c=1;
for i=contrast
   
   trialsPerContrast(c) = length(data(data(:,2)==i,2));
   
   c=c+1;
   
end
% if the number of trials for a contrast value is less than 5% then remove
% this contrast value from the data set

data_copy = data;   % I am doing it this way because otherwise inside the loop the length of data becomes smaller, Armin 3/4/2017
for i=1:length(contrast)
   
   if trialsPerContrast(i) < 0.05*length(data_copy)
      
      data(data(:,2)==contrast(i),:)=[];
      
   end
   
end



% refdefine trailsPerContrast with updated totals
contrast = unique(data(:,2))';
trialsPerContrast = nan(length(contrast),1);
c=1;
for i=contrast
   
   trialsPerContrast(c) = length(data(data(:,2)==i,2));
   
   c=c+1;
   
end


% sanity check on contrast
% take out all trials where contrast is 0
% data(data(:,2)==0,:)=[];

Data.data=data;
Data.ID=AnimalID;

%% Test weibull

% cc = unique(data(:,2))';
% nn = nan(1,length(cc));
% pp = nan(1,length(cc));
% blockID = 1;
%
% c=1;
% for istim=unique(data(:,2))'
%
%    nn(c) = length(data(data(:,2)==istim & data(:,8)==blockID,2));
%    pp(c) = length(data(data(:,2)==istim & data(:,8)==blockID & data(:,3)==1,2))/nn(c);
%
%    c=c+1;
%
% end
%
% [psychoParams, LL] = mle_fit_psycho([cc ;nn ;pp],'erf_psycho');


%% fit the model to the data
% Fval is the negative log-likelihood.

% Initialise
Nruns = length(alpha)*length(DA_val)*length(noiseSTD);
params = nan(3,1);
xanswerStore = nan(Nruns,3);
FvalStore = nan(length(alpha),length(DA_val),length(noiseSTD));
FvalStore_L = FvalStore;
FvalStore_R = FvalStore;

% Run over parameter values
counter = 0;
for k=1:length(noiseSTD)
   disp(['noiseSTD = ',num2str(k)])
   for j=1:length(DA_val)
      disp(['DA_val = ',num2str(j)])
      for i=1:length(alpha)
         params(1) = alpha(i);
         params(2) = DA_val(j);
         params(3) = noiseSTD(k);
         counter = counter + 1;
         [Fval, Fval_l, Fval_r] = FitPOMDP_GS_NLL_New(params, Data);
         xanswerStore(counter,:) = params;
         FvalStore(i,j,k) = Fval;
         FvalStore_L(i,j,k) = Fval_l;
         FvalStore_R(i,j,k) = Fval_r;
      end
   end
end


%% Individual plots

close all

set(0,'DefaultAxesPosition',[0.05,0.05,0.94,0.94])


% HISTOGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%histogram(FvalStore(:),100)


% COLOUR MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exclude extremely large values for the colour map
% redefine all Fmin values above the cut-off
% Fval_plotMax = max(FvalStore);
% FvalStoreAlt = FvalStore;
% for i=1:length(FvalStore(:))
%     
%     if FvalStore(i) > Fval_plotMax
%         
%         FvalStoreAlt(i) = Fval_plotMax;
%         
%     end
%     
% end

% Plot the colour map
% FIX ME these plots need titles and labels
PlotColourMap_NLL(alpha,DA_val,noiseSTD,FvalStore)

PlotColourMap_NLL(alpha,DA_val,noiseSTD,FvalStore_L)
PlotColourMap_NLL(alpha,DA_val,noiseSTD,FvalStore_R)
PlotColourMap_NLL(alpha,DA_val,noiseSTD,FvalStore_R+FvalStore_L)

%% 3BEST 3WORST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find the 3 best and 3 worst sets of parameters
[xanswerMin, ~] = Find3BestWorst(FvalStore,xanswerStore);

% initialise data
data_modelMin = zeros(3,length(data),22);
% data_modelMax = zeros(3,length(data),22);

% run the model with the worst and best parameters
for i=1:3
   [data_modelMin(i,:,:)] = RunPOMDP_GS_NLL(Data, xanswerMin(i,:));
   % [data_modelMax(i,:,:)] = RunPOMDP_GS_NLL(Data, xanswerMax(i,:));
end


%% Plot3BestWorst(Data,data_modelMin,data_modelMax);


% PLOT TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_model = zeros(size(data_modelMin,2),size(data_modelMin,3));
data_model(:,:) = data_modelMin(1,:,:);

% PlotTrials(Data,data_model)
PlotTrialsTogether(Data,data_model)


% LEARNING RATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotLearningRate(Data,data_model)
% PlotLearningRateReactTime(Data,data_model)


% MICE AND MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotMiceRLAndErf_GS(data,data_model)


% PLOT INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fileID = fopen([AnimalID,'_NLL.txt'],'w');
% fprintf(fileID,'%6s %6s\n','AnimalID: ',AnimalID);
% fprintf(fileID,'%6s %6s\n','Total number of trials: ',...
%    num2str(length(data)));
% fprintf(fileID,'%6s %6s\n','Number of contrast levels: ',...
%    num2str(length(trialsPerContrast)));
% fprintf(fileID,'%6s %6s\n','Mean trials per constrast: ',...
%    num2str(mean(trialsPerContrast)));
% fprintf(fileID,'%6s %6s\n','Median trials per constrast: ',...
%    num2str(median(trialsPerContrast)));
% 
% fprintf(fileID,'\n%6s\n','Model parameter values:');
% fprintf(fileID,'%6s %6s\n','alpha = ',num2str(xanswerMin(1)));
% fprintf(fileID,'%6s %6s\n','DAval = ',num2str(xanswerMin(2)));
% fprintf(fileID,'%6s %6s','noiseSTD = ',num2str(xanswerMin(3)));
% fclose(fileID);


%% ALL THE PLOTS


[~, locMin1] = min(FvalStore(:));
xanswerBest = xanswerStore(locMin1,:);

[data_modelBest] = RunPOMDP_GS_NLL(Data, [xanswerBest(1),...
   xanswerBest(2),xanswerBest(3)]);

AllThePlots(Data,data_modelBest,xanswerBest,trialsPerContrast);

