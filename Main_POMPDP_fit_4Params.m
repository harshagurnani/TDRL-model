% 1. Fit sigma, alpha, Da-val, and bias using all blocks.
% 2. If block 3/4 present, refine estimate of alpha.
% 3. Run best model on all blocks and look at psychometric and conditional
%    psychometric curves.
% 4. All Plots
% 5. Also look at performance - percent correct predicted for model (or top
%    something models).

close all

%-------------------------------------------------------------------------%
%%              Step 1: Fit parameters using all blocks
%-------------------------------------------------------------------------%

%Initialisation
% set model parameter values
alpha       = 0.03 : 0.07 : 0.7;     % learning rate
DA_val      = 0 : 0.5 : 6.0;         % dopamine value
noiseSTD    = 0.02 : 0.02 : 0.52;    % noise in belief
Qbias       = -0.2 : 0.05 : 0.2;     % bias in Q-value

% load in animal data if not already present
AnimalID='ALK045'; % 'SS031_DA' or 'SS040' or 'ALK05' or 'ALK011'
% or 'ALK028' or 'ALK017'
ExpID = '1';

if ~exist('dataAll', 'var')
  dataAll = LoadSavedDataForBehExp (AnimalID, ExpID);
end
% dataAll - untouched copy.

% copy data - no blocks removed
data =dataAll;


Data.data=data;
Data.ID=AnimalID;

%% check contrasts with insufficient number of trials (<2.5% of block trials)

contrast = unique(data(:,2))';
blocks = unique(data(:,8))';        %Should only be block 4, or none.
trialsPerContrast = nan(length(contrast)*length(unique(data(:,8))),1);
includeContrast = ones(size(trialsPerContrast));

c=1;
for b= blocks 
  for ii= contrast
    trialsPerContrast(c) = length(data(data(:,2)==ii & data(:,8)==b,2));


    if trialsPerContrast(c) < 0.05*length(data(data(:,8)==b,2))
      includeContrast(c) = 0;
    end
    c=c+1;
  end
end
%%%% HG. Not removing stimuli with lower trials - just not calculating likelihood
%%%% at those stimuli



%% find parameters sigma and bias that minimize NLL of PC
% Fval is the negative log-likelihood (NLL)

% Initialise
Nruns = length(alpha)*length(DA_val)*length(noiseSTD)*length(Qbias);
params = nan(4,1);
xanswerStore = nan(Nruns,4);
FvalStore = nan(length(alpha),length(DA_val),length(noiseSTD),length(Qbias));

% Run over parameter values
counter = 0;
iterN = 11;     %No. of iterations of model to compute probabilities
for kk=1:length(noiseSTD)
   disp(['noiseSTD = ',num2str(kk)])
   for jj=1:length(DA_val)
      disp(['DA_val = ',num2str(jj)])
      for ii=1:length(alpha)
      for bb=1:length(Qbias)
         params(1) = alpha(ii);
         params(2) = DA_val(jj);
         params(3) = noiseSTD(kk);
         params(4) = Qbias(bb);
         counter = counter + 1;
         %%%%%%%% Need to write model with bias
         [Fval, ~,~] = Get_NLL_ModelWithBias(params, Data,  includeContrast, iterN);
         xanswerStore(counter,:) = params;
         FvalStore(ii,jj,kk,bb) = Fval;
      end
      end
   end
end


% Minimize NLL - best 5 models
nModels = 5;
[xanswerMin, ~] = Find3BestWorst(FvalStore,xanswerStore, nModels);

%-------------------------------------------------------------------------%
%% Step 2: If Block 3/4 exists, fit alpha again using conditional PC
%-------------------------------------------------------------------------%
    % Fix sigma, DAval, and bias
    alpha_new       = 0.03 : 0.03 : 0.8;       % learning rate
    DA_val_new      = xanswerMin(:,2)';           % dopamine value
    noiseSTD_new    = xanswerMin(:,3)';           % noise in belief
    Qbias_new       = xanswerMin(:,4)';           % bias in QValue
    
        % Initialise
    Nruns2 = length(alpha_new);                % Only vary alpha for top 5 models (combination preserved)
    params = nan(4,1);
    xanswerStore2 = nan(Nruns2,4,nModels);
    FvalStore2 = nan(length(alpha_new),nModels);
    FvalStore2_L = FvalStore2;
    FvalStore2_R = FvalStore2;
    
if any(data(:,8)>2)
    % Block 3/4 present
    data2 = dataAll;

    % Remove other blocks
    data2( data2(:,8)==1, :) =[];
    data2( data2(:,8)==2, :) =[];
    Data2.data = data2;
    Data2.ID=AnimalID;

    % check number of trials for each contrast
    contrast2 = unique(data2(:,2))';
    trialsPerContrast2 = nan(length(contrast2)*length(unique(data2(:,8))),1);
    includeContrast2 = ones(size(trialsPerContrast2));


    c=1;
    for b=unique(data2(:,8))'   
      for ii=contrast2
        trialsPerContrast2(c) = length(data2(data2(:,2)==ii & data2(:,8)==b,2));


        if trialsPerContrast2(c) < 0.025*length(data2(data2(:,8)==b,2))
        includeContrast2(c) = 0;
        end
        c=c+1;
      end
    end
    




    iterN = 21;     
    % Run over parameter values
    
    for kk=1:nModels
          counter = 0; % Restart counting number of runs for each model
          % No. of runs --> ii
          for ii=1:length(alpha_new)
             params(1) = alpha_new(ii);
             % All others taken from top models
             params(2) = DA_val_new(kk);
             params(3) = noiseSTD_new(kk);
             params(4) = Qbias_new(kk);
             counter = counter + 1;
             [Fval, Fval_l, Fval_r] = Get_NLL_ModelWithBias(params, Data2, includeContrast2, iterN);
             xanswerStore2(counter,:,kk) = params;
             FvalStore2(ii,kk) = Fval;
             FvalStore2_L(ii,kk) = Fval_l;
             FvalStore2_R(ii,kk) = Fval_r;
          end
    end
    
    %% Find the best 5 alpha (+other parameters)
    [xanswerMin2, ~] = Find3BestWorst(FvalStore2_L+FvalStore2_R,xanswerStore2,1,true );
    xanswerMin2 = squeeze(permute( xanswerMin2, [3 2 1]));      %1 alpha for each of top 5 models --> 5 new top models
    refinedAlpha = 1;
else
    xanswerMin2 = xanswerMin;
    refinedAlpha = 0;
end
% 
% % %% Plot psychometric curves of fitted data
% % PlotColourMap_NLL(alpha_new,DA_val_new,noiseSTD_new,FvalStore2, 'NLL of Block-wise psychometric curves')
% % PlotColourMap_NLL(alpha_new,DA_val_new,noiseSTD_new,FvalStore2_L+FvalStore2_R, 'NLL of Conditional P(R)')
% 
% 

%-------------------------------------------------------------------------%
%% Step 4: Run model on all blocks to exaine Psychometric curves
%-------------------------------------------------------------------------%
% 
data3     =dataAll;
Data3.data=data3;
Data3.ID  =AnimalID;

data_models = nan(length(data3), 22, nModels);
action      = nan(length(data3), 99, nModels);
correct     = action;

% %%%% HG. Not removing stimuli with lower trials - just not calculating likelihood
% %%%% at those stimuli
% 

for model = 1:nModels
    params = xanswerMin2(model,:);
% xanswerBest2 = xanswerMin2(1,:);
% 
% %% Plot all psychometric and conditional psychometric curves
[data_models(:,:,model), action(:,:,model), correct(:,:,model), ~, ~ ] = RunPOMDP_GS_NLL_returnDetails(Data3,params);
end
% PlotPC_Mice_Model(Data2.data, data_modelBest2, action2, includeContrast2)
% PlotCondPC_Mice_Model(Data2.data, data_modelBest2, action2, correct2, includeContrast2)
% 
% 
% 
% %-------------------------------------------------------------------------%
% %% Step 5: Plot percent predicted [as previously only P(R|correct) was fit.]
% %-------------------------------------------------------------------------%