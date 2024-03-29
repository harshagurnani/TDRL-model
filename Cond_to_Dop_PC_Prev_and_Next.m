function [ stimuli, per_mice_next, per_mice_prev ] = Cond_to_Dop_PC_Prev_and_Next(data_mice, varargin)

% changing y-axis values
if any(data_mice(:,3) == -1)
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
end
contrasts = unique(data_mice(:,2))';


if nargin>2
    blocks = varargin{2};
else
    blocks=unique(data_mice(:,8))';
end

if nargin >1
    includeContrast = varargin{1}; 
else
    includeContrast = ones(length(contrasts),length(blocks));
end

includeContrast = reshape(includeContrast, length(contrasts), length(blocks) );


per_mice_next = cell( length(blocks),1); per_mice_prev = per_mice_next;
stimuli = cell( length(blocks),1);

%%% need to change averaging across both/none blocks because P(R) is not
%%% the right measure.
b=0;
DopTrial = nan(length(data_mice),1);
for blockID = blocks
   b=b+1;
   if blockID == 1
       DChoice = 0;
   elseif blockID == 2
       DChoice = 1;
   elseif blockID == 3
       DChoice = [0 1];
   elseif blockID == 4
       DChoice = 2;
   end
   
   
   for jj=1:length(data_mice)
      DopTrial(jj,1) = any(data_mice(jj,3)==DChoice); 
   end
   
   per_mice_next{b} = nan(2,length(contrasts));
   % Top row is correct, bottom row is error.
   per_mice_prev{b} = per_mice_next{b};

   c=1;
   for istim = contrasts
      
      %After/Before correct dopamine-reward choices.
      % Each stimulus/block plus previous action was correct +
      % dopamine-rewarding
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & data_mice(2:end,7) < 5 & ...
          DopTrial(1:end-1,1) & data_mice(1:end-1,10)==1 ];
      per_mice_next{b}(1,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & data_mice(1:end-1,7) < 5 & ...
          DopTrial(2:end,1) & data_mice(2:end,10)==1 ; false];
      per_mice_prev{b}(1,c)= nanmean(data_mice(id,3)) ;
      
      %After/Before error trials where the action is associated with dopamine on correct trials.
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & data_mice(2:end,7) < 5 & ...
          DopTrial(1:end-1,1) & data_mice(1:end-1,10)==0 ];
      per_mice_next{b}(2,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & data_mice(1:end-1,7) < 5 & ...
          DopTrial(2:end,1) & data_mice(2:end,10)==0 ; false];
      per_mice_prev{b}(2,c)= nanmean(data_mice(id,3)) ;
      
      c=c+1;
   end
   
   ic = ( includeContrast(:,b)==1 );
   stimuli{b} = contrasts(ic); 
   per_mice_next{b} = per_mice_next{b}(:,ic);
   per_mice_prev{b} = per_mice_prev{b}(:,ic);
   
end


end

