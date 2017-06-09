function [ stimuli, shifts, per_mice_next, includeContrast ] = SingleTrialShift_EasyDiffStim_with_Dop( dataGiven , varargin )
% previous Stimulus is easy (highest contrast) or difficult (lowest contrast)
% After correct actions that give dopamine
% return actual PC and includeContrast if needed

if nargin == 3
    % Multiple iters of model
      action = varargin{1};
      correct = varargin{2};
      nIters = size(action,2);
    else
      nIters = 1;
      action = dataGiven(:,3);
      correct = dataGiven(:,10);
end
    

% changing choices to 0 and 1
if any(dataGiven(:,3) == -1)
  dataGiven(:,3) = (1 + dataGiven(:,3)) ./ 2;
end
if any(action == -1)
    action = (1+action)/2;
end


contrasts = unique(dataGiven(:,2))';


blocks= [1 2];%unique(data(:,8))';

% Which stimuli to include (at least 5% frequency)
includeContrast = ones(length(contrasts),1);
trialsPerContrast = nan(length(contrasts), 1);
c=1;
for ii=contrasts
    trialsPerContrast(c) = length(dataGiven(dataGiven(:,2)==ii,2));
    if trialsPerContrast(c) < 0.05*length(dataGiven)
    includeContrast(c) = 0;
    end
    c=c+1;
end


% Contrasts, DA-L or DA-R, Easy or Diff


stimuli = contrasts(includeContrast ==1 );
ic = includeContrast==1;

easystim = [ min(stimuli)                max(stimuli) ];           %Highest contrasts
diffstim = [ max(stimuli(stimuli<0))  0  min(stimuli(stimuli>0)) ];%Lowest contrasts

IsStimEasy = (dataGiven(:,2)<=easystim(1) | dataGiven(:,2)>=easystim(2));
IsStimDiff = (dataGiven(:,2)>=1.2*diffstim(1) & dataGiven(:,2)<=diffstim(2));

%%% checking shape of array
IsStimEasy = reshape( IsStimEasy, [ length(dataGiven), 1 ]);
IsStimDiff = reshape( IsStimDiff, [ length(dataGiven), 1 ]);


per_mice_next = nan(length(contrasts), 2, 2, nIters); 
per_mice_prev = per_mice_next;
shifts = nan(length(stimuli),       2,               2,           nIters);
%           current stimulus, block number,   previous stimulus,  nIters
% Dim 3:    Indx 1 is for easy stim, indx 2 is for diff stim.

data = dataGiven;       % Copy used and changed
for iter = 1:nIters
    
    
    data(:,3) = action(:,iter);
    data(:,10) = correct(:, iter);
    
    
    b=0;
    for blockID = blocks
       b=b+1;

       if blockID == 1
           DChoice = 0;
       elseif blockID == 2
           DChoice = 1;
       end
       DopTrial =  (data(:,3) == DChoice);
       DopTrial = reshape( DopTrial, [ length(dataGiven), 1 ]);
       
       c=1;
       for istim = contrasts

          %After/Before easy stim + dopamine.

          id= [false;data(2:end,2)==istim  & data(2:end,8)==blockID & data(2:end,7) < 5 & ... current trial
              DopTrial(1:end-1,1) & data(1:end-1,10)==1  & IsStimEasy(1:end-1)];            % conditioning trial
           per_mice_next(c,blockID,1,iter)= nanmean(data(id,3)) ;
          id= [data(1:end-1,2)==istim & data(1:end-1,8)==blockID & data(1:end-1,7) < 5 & ...
             DopTrial(2:end,1) & data(2:end,10)==1       & IsStimEasy(2:end); false];
           per_mice_prev(c,blockID,1,iter)= nanmean(data(id,3)) ;

          %After/Before difficult stim + dopamine
          id= [false;data(2:end,2)==istim & data(2:end,8)==blockID & data(2:end,7) < 5 & ...
              DopTrial(1:end-1,1) & data(1:end-1,10)==1  & IsStimDiff(1:end-1)];
           per_mice_next(c,blockID,2)= nanmean(data(id,3)) ;
          id= [data(1:end-1,2)==istim & data(1:end-1,8)==blockID & data(1:end-1,7) < 5 & ...
              DopTrial(2:end,1) & data(2:end,10)==1      & IsStimDiff(2:end); false];
           per_mice_prev(c,blockID,2,iter)= nanmean(data(id,3)) ;

          c=c+1;
       end


       shifts(:, blockID, :, iter) = per_mice_next(ic,blockID,:,iter)- per_mice_prev(ic,blockID,:,iter);

    end

end

end

