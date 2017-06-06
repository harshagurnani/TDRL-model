function [ stimuli, shifts ] = SingleTrialShift_StimLevel_and_Dop( data_mice , varargin )

markerShape = {'o', 's', '^', 'v'};
blockColor = {'b','r','g','k'};
blockLabel = {'DA-L', 'DA-R', 'DA-Both', 'DA-None'};

% changing y-axis values
if any(data_mice(:,3) == -1)
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
end
contrasts = unique(data_mice(:,2))';


if nargin>2
    blocks = varargin{2};
else
    blocks= [1 2];%unique(data_mice(:,8))';
end

if nargin >1
    includeContrast = varargin{1}; 
else
    includeContrast = ones(length(contrasts),length(blocks));
    trialsPerContrast = nan(length(contrasts), length(blocks));

    c=1;b=0;
    for blockID=blocks  
        b=b+1;
      for ii=contrasts
        trialsPerContrast(c) = length(data_mice(data_mice(:,2)==ii & data_mice(:,8)==blockID,2));


        if trialsPerContrast(c) < 0.05*length(data_mice(data_mice(:,8)==blockID,2))
        includeContrast(c) = 0;
        end
        c=c+1;
      end
    end
    

end

includeContrast = reshape(includeContrast, length(contrasts), length(blocks) );

per_mice_next = cell( length(blocks),1); per_mice_prev = per_mice_next;
shifts = cell( length(blocks),1);
stimuli = cell( length(blocks),1);

figure; hold on;
h = zeros(2*length(blocks),1);
mLabel = cell(2*length(blocks),1);

easy_stim = nan(length(data_mice),1);

b=0;
for blockID = blocks
   b=b+1;
   stims = contrasts( includeContrast(:,b)==1);
   stimlevels = sort( stims, 'ascend');
   diffStim = [ stimlevels(2) stimlevels(end-1) ]
%    diffStim = 
2*[max(stims(stims<0)) min(stims(stims>0))]
   if blockID == 1
       DChoice = 0;
   elseif blockID == 2
       DChoice = 1;
   elseif blockID == 3
       DChoice = [0 1];
   end
   
   per_mice_next{b} = nan(2,length(contrasts));
   % Top row is after easy, bottom row is after diff stim.
   per_mice_prev{b} = per_mice_next{b};

   c=1;
   for istim = contrasts

      %After/Before easy stim+dopamine.
      for tt=1:length(data_mice)
         easy_stim(tt) = (data_mice(tt,2)<diffStim(1) || data_mice(tt,2)>diffStim(2));
      end
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          any(data_mice(1:end-1,3)==DChoice) & easy_stim(1:end-1) & ...
          data_mice(1:end-1,10)==1 & data_mice(2:end,7) < 5];
      per_mice_next{b}(1,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & ...
          any(data_mice(2:end,3)==DChoice) & easy_stim(2:end) ...
          & data_mice(2:end,10)==1 & data_mice(1:end-1,7) < 5; false];
      per_mice_prev{b}(1,c)= nanmean(data_mice(id,3)) ;
      
      %After/Before difficult stim + dopamine
      id=[false;data_mice(2:end,2)==istim & data_mice(2:end,8)==blockID & ...
          any(data_mice(1:end-1,3)==DChoice) & ~easy_stim(1:end-1) ...
          & data_mice(1:end-1,10)==1 & data_mice(2:end,7) < 5];
      per_mice_next{b}(2,c)= nanmean(data_mice(id,3)) ;
      id=[data_mice(1:end-1,2)==istim & data_mice(1:end-1,8)==blockID & ...
          any(data_mice(2:end,3)==DChoice) & ~easy_stim(2:end) ...
          & data_mice(2:end,10)==1 & data_mice(1:end-1,7) < 5; false];
      per_mice_prev{b}(2,c)= nanmean(data_mice(id,3)) ;
      
      c=c+1;
   end
   
   ic = ( includeContrast(:,b)==1 );
   stimuli{b} = contrasts(ic); 
   per_mice_next{b} = per_mice_next{b}(:,ic);
   per_mice_prev{b} = per_mice_prev{b}(:,ic);
   shifts{b} = per_mice_next{b} - per_mice_prev{b};
   
   h(2*b-1) = plot( stimuli{b}, per_mice_next{b}(1,:), 'marker', markerShape{blockID},...
                    'color', blockColor{blockID}, 'linestyle','--');
   
   h(2*b) = plot( stimuli{b}, per_mice_next{b}(2,:), 'marker', markerShape{blockID},...
                    'color', blockColor{blockID}, 'linestyle','-','markerfacecolor', blockColor{blockID});
   mLabel{2*b-1} = sprintf('Block %s - After Easy Stimulus', blockLabel{blockID});                          
   mLabel{2*b} = sprintf('Block %s - After Diff Stimulus', blockLabel{blockID});                          
end

legend(h, mLabel{:})


end

