function [ stimuli, shifts ] = SingleTrialShift_StimLevel_and_Correct( data_mice , varargin )

markerShape = {'o', 's', '^', 'v'};
blockColor = {'b','r','g','k'};
blockLabel = {'Left', 'Right', 'DA-Both', 'DA-None'};

% changing y-axis values
if any(data_mice(:,3) == -1)
  data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
end
contrasts = unique(data_mice(:,2))';



blocks= [1 2];%unique(data_mice(:,8))';



includeContrast = ones(length(contrasts),1);
trialsPerContrast = nan(length(contrasts), 1);

c=1;

for ii=contrasts
    trialsPerContrast(c) = length(data_mice(data_mice(:,2)==ii,2));
    if trialsPerContrast(c) < 0.05*length(data_mice)
    includeContrast(c) = 0;
    end
    c=c+1;
end

per_mice_next = nan(length(contrasts),2,2); 
per_mice_prev = per_mice_next;
% Contrasts, DA-L or DA-R, Easy or Diff


stimuli = contrasts(includeContrast ==1 );
ic = includeContrast==1;
figure; hold on;
h = zeros(4,1);
mLabel = cell(4,1);

IsStimEasy = nan(length(data_mice),1);
IsStimDiff = IsStimEasy;
stimSorted = sort( stimuli, 'ascend');

easystim = [ stimSorted(1) stimSorted(end) ];
diffStim = [ max(stimuli(stimuli<0)) 0 min(stimuli(stimuli>0)) ];

IsStimEasy = (data_mice(:,2)<=diffStim(1) | data_mice(:,2)>=diffStim(2));

shifts = nan(length(stimuli), 2,2);

b=0;
for blockID = blocks
   b=b+1;
   
   if blockID == 1
       DChoice = 0;
   elseif blockID == 2
       DChoice = 1;
   end
   
  
   % Indx 1 is after easy,indx 2 is after diff stim.

   c=1;
   for istim = contrasts

      %After/Before easy stim+dopamine.
      
      id=[false;data_mice(2:end,2)==istim  & ...
          data_mice(1:end-1,3)==DChoice & IsStimEasy(1:end-1) & ...
          data_mice(1:end-1,10)==1 & data_mice(2:end,7) < 5];
      per_mice_next(c,blockID,1)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & ...
          data_mice(2:end,3)==DChoice & IsStimEasy(2:end) ...
          & data_mice(2:end,10)==1 & data_mice(1:end-1,7) < 5; false];
      per_mice_prev(c,blockID,1)= nanmean(data_mice(id,3)) ;
      
      %After/Before difficult stim + dopamine
      id=[false;data_mice(2:end,2)==istim & ...
          data_mice(1:end-1,3)==DChoice & ~IsStimEasy(1:end-1) ...
          & data_mice(1:end-1,10)==1 & data_mice(2:end,7) < 5];
      per_mice_next(c,blockID,2)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & ...
          data_mice(2:end,3)==DChoice & ~IsStimEasy(2:end) ...
          & data_mice(2:end,10)==1 & data_mice(1:end-1,7) < 5; false];
      per_mice_prev(c,blockID,2)= nanmean(data_mice(id,3)) ;
      
      c=c+1;
   end
   
   
   shifts(:, blockID, :) = per_mice_next(ic,blockID,:)- per_mice_prev(ic,blockID,:);
   
%    shifts{b} = per_mice_next{b} - per_mice_prev{b};
   
   h(2*b-1) = plot( stimuli, per_mice_next(ic,blockID,1), 'marker', markerShape{blockID},...
                    'color', blockColor{blockID}, 'linestyle','--');
   
   h(2*b) = plot( stimuli, per_mice_next(ic,blockID,2), 'marker', markerShape{blockID},...
                    'color', blockColor{blockID}, 'linestyle','-','markerfacecolor', blockColor{blockID});
   mLabel{2*b-1} = sprintf(' %s - After Easy Stimulus', blockLabel{blockID});                          
   mLabel{2*b}   = sprintf(' %s - After Diff Stimulus', blockLabel{blockID});                          
end

legend(h, mLabel{:})


end

