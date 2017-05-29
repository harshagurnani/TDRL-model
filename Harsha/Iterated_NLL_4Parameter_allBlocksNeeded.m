% 1. Fit sigma and bias using block 4/3 by maximising likelihood of psychometric curve.
% 2. Fix sigma and bias, and run again on same data to fit alpha by
%    maximising likelihood of conditional psychometric curve ==> Do it on
%    Block 4 only?
% 3. Fit Dopamine value using Block 3? or fix some value./ Minimize NLL of
%    conditional PC
% 4. Run model on all blocks and look at psychometric and conditional
%    psychometric curves.
% 5. Also look at performance - percent correct predicted for model (or top
%    something models).

close all

%-------------------------------------------------------------------------%
%%              Step 1: Fit sigma and bias
%-------------------------------------------------------------------------%

%Initialisation
% set model parameter values
alpha       = 0.0 : 0.1 : 0.7;       % learning rate
DA_val      = 0.0;% : 0.5 : 6.0;     % dopamine value
noiseSTD    = 0.02 : 0.02 : 0.42;    % noise in belief
Qbias       = -0.2 : 0.05 : 0.2;     %bias in Q-value

% load in animal data if not already present
AnimalID='SS037'; % 'SS031_DA' or 'SS040' or 'ALK05' or 'ALK011'
% or 'ALK028' or 'ALK017'
ExpID = '1';

if ~exist('dataAll', 'var')
  dataAll = LoadSavedDataForBehExp (AnimalID, ExpID);
end

% copy data and remove all blocks except No DA.
data =dataAll;
data(data(:,8)==3,:)=[];
data(data(:,8)==1,:)=[];
data(data(:,8)==2,:)=[];

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


    if trialsPerContrast(c) < 0.025*length(data(data(:,8)==b,2))
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
iterN = 21;     %No. of iterations of model to compute probabilities
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


% Fix sigma and bias by minimizing NLL
[xanswerMin, ~] = Find3BestWorst(FvalStore,xanswerStore,1);

%-------------------------------------------------------------------------%
%% Step 2: Run again to get conditional psychometric curves fitted
%-------------------------------------------------------------------------%

% Fix noise and sigma
alpha_new       = 0.03 : 0.03 : 1;         % learning rate
DA_val_new      = 0;                       % dopamine value
noiseSTD_new    = xanswerMin(3);           % noise in belief
Qbias_new       = xanswerMin(4);           % bias in QValue

% Initialise
Nruns2 = length(alpha_new)*length(DA_val_new)*length(noiseSTD_new)*length(Qbias_new);
params = nan(4,1);
xanswerStore2 = nan(Nruns2,4);
FvalStore2 = nan(length(alpha_new), length(DA_val_new), length(noiseSTD_new), length(Qbias_new));
FvalStore2_L = FvalStore2;
FvalStore2_R = FvalStore2;

iterN = 21;     
% Run over parameter values
counter = 0;
for kk=1:length(noiseSTD_new)
   disp(['noiseSTD = ',num2str(kk)])
   for jj=1:length(DA_val_new)
      disp(['DA_val = ',num2str(jj)])
      for ii=1:length(alpha_new)
      for bb=1:length(Qbias_new)
         params(1) = alpha_new(ii);
         params(2) = DA_val_new(jj);
         params(3) = noiseSTD_new(kk);
         params(4) = Qbias_new(bb);
         counter = counter + 1;
         [Fval, Fval_l, Fval_r] = Get_NLL_ModelWithBias(params, Data, includeContrast, iterN);
         xanswerStore2(counter,:) = params;
         FvalStore2(ii,jj,kk,bb) = Fval;
         FvalStore2_L(ii,jj,kk,bb) = Fval_l;
         FvalStore2_R(ii,jj,kk,bb) = Fval_r;
      end
      end
   end
end

%% Plot psychometric curves of fitted data
PlotColourMap_NLL(alpha_new,DA_val_new,noiseSTD_new,FvalStore2, 'NLL of Block-wise psychometric curves')
PlotColourMap_NLL(alpha_new,DA_val_new,noiseSTD_new,FvalStore2_L+FvalStore2_R, 'NLL of Conditional P(R)')

%% Find the 3 best sets of parameters
[xanswerMin2, ~] = Find3BestWorst(FvalStore2_L+FvalStore2_R,xanswerStore2);

%-------------------------------------------------------------------------%
%% Step 3: Fit dopamine value using block 3. If absent, just set manually
%-------------------------------------------------------------------------%
xanswerBest = xanswerMin2(1,:);

data3 =dataAll;
Data3.data=data3;
Data3.ID=AnimalID;
data3(data3(:,8)==1,:)=[];
data3(data3(:,8)==2,:)=[];
data3(data3(:,8)==4,:)=[];

if any(data3(:,8) == 3)
%% If Block 3 present, fit DA using that

    alpha3       = xanswerBest(1);       % learning rate
    DA_val3      = 0.0 : 0.5 : 6.0;      % dopamine value
    noiseSTD3    = xanswerBest(3);       % noise in belief
    Qbias3       = xanswerBest(4);       % bias in Q-value
    
    
    % check no. of trials per contrast
    contrast = unique(data3(:,2))';
    blocks   = 3;
    trialsPerContrast3 = nan(length(contrast),1);
    includeContrast3 = ones(size(trialsPerContrast3));

    c=1;
    for b= blocks 
      for ii= contrast
        trialsPerContrast3(c) = length(data3(data3(:,2)==ii & data3(:,8)==b,2));


        if trialsPerContrast3(c) < 0.025*length(data3(data3(:,8)==b,2))
          includeContrast3(c) = 0;
        end
        c=c+1;
      end
    end
    
    % Initialise
    Nruns3 = length(alpha3)*length(DA_val3)*length(noiseSTD3)*length(Qbias3);
    params = nan(4,1);
    xanswerStore3 = nan(Nruns3,4);
    FvalStore3 = nan(length(alpha3),length(DA_val3),length(noiseSTD3),length(Qbias3));
    FvalStore3_L = FvalStore3;
    FvalStore3_R = FvalStore3;
    
    % Run over parameter values
    counter = 0;
    iterN = 21;     %No. of iterations of model to compute probabilities
    for kk=1:length(noiseSTD3)
       disp(['noiseSTD = ',num2str(kk)])
       for jj=1:length(DA_val3)
          disp(['DA_val = ',num2str(jj)])
          for ii=1:length(alpha3)
          for bb=1:length(Qbias3)
             params(1) = alpha3(ii);
             params(2) = DA_val3(jj);
             params(3) = noiseSTD3(kk);
             params(4) = Qbias3(bb);
             counter = counter + 1;
             %%%%%%%% Need to write model with bias
             [Fval, Fval_L, Fval_R] = Get_NLL_ModelWithBias(params, Data3,  includeContrast3, iterN);
             xanswerStore3(counter,:) = params;
             FvalStore3(ii,jj,kk,bb) = Fval;
             FvalStore3_L(ii,jj,kk,bb) = Fval_L;
             FvalStore3_R(ii,jj,kk,bb) = Fval_R;
          end
          end
       end
    end


    % Fix DA-val (everything else already fixed)
    [xanswerMin3, ~] = Find3BestWorst(FvalStore3_L+FvalStore3_R,xanswerStore3,1);
    xanswerBest = xanswerMin3;
    
else
    xanswerBest(2) = 3;

end

%-------------------------------------------------------------------------%
%% Step 4: Run model on all blocks to exaine Psychometric curves
%-------------------------------------------------------------------------%

data2 =dataAll;
Data2.data=data2;
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
%%%% HG. Not removing stimuli with lower trials - just not calculating likelihood
%%%% at those stimuli

xanswerBest2 = xanswerBest;

%% Plot all psychometric and conditional psychometric curves
[data_modelBest2, action2, correct2, QL2, QR2 ] = RunPOMDP_GS_NLL_returnDetails(Data2, xanswerBest2);

PlotPC_Mice_Model(Data2.data, data_modelBest2, action2, includeContrast2)
PlotCondPC_Mice_Model(Data2.data, data_modelBest2, action2, correct2, includeContrast2)



%-------------------------------------------------------------------------%
%% Step 5: Plot percent predicted [as previously only P(R|correct) was fit.]
%-------------------------------------------------------------------------%