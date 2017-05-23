close all

% TESTING2 set model parameter values
alpha       = 0.0 : 0.1 : 1;       % learning rate
DA_val      = 0.0 : 0.5 : 6.0;       % dopamine value
noiseSTD    = 0.02 : 0.02 : 0.42;    % noise in belief

% load in animal data
AnimalID='ALK028'; % 'SS031_DA' or 'SS040' or 'ALK05' or 'ALK011'
% or 'ALK028' or 'ALK017'

ExpID = '1';

if ~exist(dataAll, 'var')
dataAll = LoadSavedDataForBehExp (AnimalID, ExpID);
end

% % blockIDs 1 & 2 represent asymmetric dopamine reward, not included for
% estimating sigma/alpha
% data(data(:,8)==3,:)=[];
% data(data(:,8)==4,:)=[];

data(data(:,8)==1,:)=[];
data(data(:,8)==2,:)=[];


% check number of trials for each contrast
contrast = unique(data(:,2))';
trialsPerContrast = nan(length(contrast),1);
c=1;
for i=contrast
   
   trialsPerContrast(c) = length(data(data(:,2)==i,2));
   
   c=c+1;
   
end

%%%% HG. Not removing stimuli with lower trials - just not calculating likelihood
%%%% at those stimuli

Data.data=data;
Data.ID=AnimalID;

%% fit the model to the data
% Fval is the negative log-likelihood.

% Initialise
Nruns = length(alpha)*length(DA_val)*length(noiseSTD);
params = nan(3,1);
xanswerStore = nan(Nruns,3);
FvalStore = nan(length(alpha),length(DA_val),length(noiseSTD));

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
         %%%%%%%% Need to exclude trials for this fit
         [Fval, ~,~] = FitPOMDP_GS_NLL_New(params, Data);
         xanswerStore(counter,:) = params;
         FvalStore(i,j,k) = Fval;
      end
   end
end


%% Individual plots for old data
set(0,'DefaultAxesPosition',[0.05,0.05,0.94,0.94])
PlotColourMap_NLL(alpha,DA_val,noiseSTD,FvalStore)

%% 3BEST for limiting DA_val and noise_STD

[xanswerMin, ~] = Find3BestWorst(FvalStore,xanswerStore,5);

% %% Run again to get conditional psychometric curves fitted
% 
% % Fix noise
% alpha_new       = 0.03 : 0.03 : 1.2;        % learning rate
% DA_val_new      = sort(xanswerMin(:,2));    % dopamine value
% noiseSTD_new    = xanswerMin(1,3);          % noise in belief
% 
% % Initialise
% Nruns2 = length(alpha_new)*length(DA_val_new)*length(noiseSTD_new);
% params = nan(3,1);
% xanswerStore2 = nan(Nruns2,3);
% FvalStore2 = nan(length(alpha_new),length(DA_val_new),length(noiseSTD_new));
% FvalStore2_L = FvalStore2;
% FvalStore2_R = FvalStore2;
% 
% iterN = 50;     %No. of iterations
% % Run over parameter values
% counter = 0;
% for k=1:length(noiseSTD_new)
%    disp(['noiseSTD = ',num2str(k)])
%    for j=1:length(DA_val_new)
%       disp(['DA_val = ',num2str(j)])
%       for i=1:length(alpha_new)
%          params(1) = alpha_new(i);
%          params(2) = DA_val_new(j);
%          params(3) = noiseSTD_new(k);
%          counter = counter + 1;
%          [Fval, Fval_l, Fval_r] = FitPOMDP_GS_NLL_New(params, Data, iterN);
%          xanswerStore2(counter,:) = params;
%          FvalStore2(i,j,k) = Fval;
%          FvalStore2_L(i,j,k) = Fval_l;
%          FvalStore2_R(i,j,k) = Fval_r;
%       end
%    end
% end
% 
% %% New plots
% 
% PlotColourMap_NLL(alpha_new,DA_val_new,noiseSTD_new,FvalStore2, 'NLL of Block-wise psychometric curves')
% PlotColourMap_NLL(alpha_new,DA_val_new,noiseSTD_new,FvalStore2_L, 'NLL of Fraction Right-Choices after Correct Left choices')
% PlotColourMap_NLL(alpha_new,DA_val_new,noiseSTD_new,FvalStore2_R, 'NLL of Fraction Right-Choices after Correct Right choices')
% 
% %% 3BEST 
% 
% % find the 3 best sets of parameters
% [xanswerMin2, ~] = Find3BestWorst(FvalStore2_L+FvalStore2_R,xanswerStore2);
% 
% %% Run model for fitted parameters
% 
% % initialise data
% data_modelMin = zeros(3,length(data),22);
% % data_modelMax = zeros(3,length(data),22);
% 
% % run the model with the worst and best parameters
% for i=1:3
%    [data_modelMin(i,:,:),~,~] = RunPOMDP_GS_NLL(Data, xanswerMin2(i,:));
%    % [data_modelMax(i,:,:)] = RunPOMDP_GS_NLL(Data, xanswerMax(i,:));
% end
% 
% 
% %% Plot3BestWorst(Data,data_modelMin,data_modelMax);
% 
% 
% % PLOT TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% data_model = zeros(size(data_modelMin,2),size(data_modelMin,3));
% 
% xanswerBest = xanswerMin2(1,:);
% 
% [data_modelBest, action, correct, QL,QR] = RunPOMDP_GS_NLL_returnDetails(Data, [xanswerBest(1),...
%    xanswerBest(2),xanswerBest(3)]);
% 
% %%
% AllThePlots(Data,data_modelBest,xanswerBest,trialsPerContrast);
% 
% %% Cond PC
% PlotCondPC_Mice_Model(Data.data, data_modelBest, action, correct)
% 
% %% Plot Q-Values
% PlotQValues(data_modelBest, action, correct, QR, QL)
