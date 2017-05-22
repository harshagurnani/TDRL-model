alpha =     0.4;
DA_val =    0;
noiseSTD = 0.3;

shiftBF = [0 0.05 0.1 0.15 0.2 0.3 0.4];
noiseL = noiseSTD*[0.5 0.75 1 1.2 1.5 1.8 2];
noiseR = noiseSTD*1;

AnimalID='ALK022'; 
ExpID = '1';

% data = LoadSavedDataForBehExp(AnimalID, ExpID);

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
data_model = zeros(14,length(data),22);



% run the model
nparams=7;
for jj=1:nparams
    
 params_shift = [alpha, DA_val, noiseSTD, shiftBF(jj)]; 
 params_noisy = [alpha, DA_val, noiseSTD, noiseL(jj), noiseR ];
 data_model(jj,:,:) = RunPOMDP_GS_NLL_shiftBelief(Data, params_shift);
 data_model(nparams+jj,:,:) = RunPOMDP_GS_NLL_shiftBelief(Data, params_noisy);
 
end

data1= data; data2=data;
data1(:,3) = data_model(1,:,22);
data2(:,3) = data_model(10,:,22);
 %% Plotting
 
 for jj=1:nparams
     PlotMiceRLAndErf_GS(data1,squeeze(data_model(jj,:,:)))
 end
 
 
 for jj=1:nparams
     PlotMiceRLAndErf_GS(data2,squeeze(data_model(nparams+jj,:,:)))
 end
 
