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

b=0;
for blockID = blocks
   b=b+1;
   if blockID == 1
       DChoice = 0;
   elseif blockID == 2
       DChoice = 1;
   elseif blockID == 3
       DChoice = [0 1];
   end
   
   per_mice_next{b} = nan(2,length(contrasts));
   % Top row is correct, bottom row is error.
   per_mice_prev{b} = per_mice_next{b};

   c=1;
   for istim = contrasts
      
      %After/Before correct dopamine-reward choices.
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          any(data_mice(1:end-1,3)==DChoice) & data_mice(1:end-1,10)==1 & data_mice(2:end,7) < 5];
      per_mice_next{b}(1,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & ...
          any(data_mice(1:end-1,3)==DChoice) & data_mice(2:end,10)==1 & data_mice(1:end-1,7) < 5; false];
      per_mice_prev{b}(1,c)= nanmean(data_mice(id,3)) ;
      
      %After/Before error trials where the action is associated with dopamine on correct trials.
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          any(data_mice(1:end-1,3)==DChoice) & data_mice(1:end-1,10)==0 & data_mice(2:end,7) < 5];
      per_mice_next{b}(2,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & ...
          any(data_mice(1:end-1,3)==DChoice) & data_mice(2:end,10)==0 & data_mice(1:end-1,7) < 5; false];
      per_mice_prev{b}(2,c)= nanmean(data_mice(id,3)) ;
      
      c=c+1;
   end
   
   ic = ( includeContrast(:,b)==1 );
   stimuli{b} = contrasts(ic); 
   per_mice_next{b} = per_mice_next{b}(:,ic);
   per_mice_prev{b} = per_mice_prev{b}(:,ic);
   
end


end

