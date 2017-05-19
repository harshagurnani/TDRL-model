alpha =     0.4;
DA_val =    2;
noiseSTD = 0.3;

shiftBF = 0.25;
noiseL = noiseSTD*2;
noiseR = noiseSTD*1;

AnimalID='ALK027'; 
ExpID = '22';

data = LoadSavedDataForBehExp(AnimalID, ExpID);

data(data(:,8)==3,:)=[];
data(data(:,8)==4,:)=[];

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

%% Run model for fitted parameters

% initialise data
data_modelMin = zeros(2,length(data),22);


params_shift = [alpha, DA_val, noiseSTD, shiftBF];
params_noisy = [alpha, DA_val, noiseSTD, noiseL, noiseR ];

% run the model
 data_modelMin(1,:,:) = RunPOMDP_GS_NLL_shiftBelief(Data, params_shift);
 data_modelMin(2,:,:) = RunPOMDP_GS_NLL_shiftBelief(Data, params_noisy);
 
 %% Plotting
 
 

 
PlotMiceRLAndErf_GS(data,squeeze(data_modelMin(1,:,:)))
PlotMiceRLAndErf_GS(data,squeeze(data_modelMin(2,:,:)))
