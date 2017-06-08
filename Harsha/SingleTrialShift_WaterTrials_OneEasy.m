function [ shifts ] = SingleTrialShift_WaterTrials_OneEasy( dataGiven , varargin )

    if nargin == 3
    % Multiple iters of model
      action = varargin{1};
      correct = varargin{2};
      numRuns = size(action,2);
    else
      numRuns = 1;
      action = dataGiven(:,3);
      correct = dataGiven(:,10);
    end
    blocks = [1 2];
    contrast = unique(dataGiven(:,2))';
    
    % Contrasts to include (sufficient number of trials)
    trialsPerContrast = nan(length(contrast),1);
    

    
    
    % 2 = Easy , 1 = Difficult (second lowest stimulus), 0 = zero contrast
    stimuli = [ -2  -1  0  1  2 ];
    shifts = nan( length(stimuli), 2, length(blocks), numRuns );
    shift_n = nan(length(stimuli),2, length(blocks), numRuns);
    % stimulus, conditional to left/right actions, block, animalID
    for AID=1:numRuns
        dataAll = dataGiven;
        % Replace actions/correct with iteration values, if applicable
        dataAll(:,3) = action(:,AID);
        dataAll(:,10)   = correct(:,AID);
        includeContrast = ones(length(contrast),1);
        b=0;
        for blockID = blocks
        b=b+1;
        %% Get contrats to include
        data=dataAll; 
        data( data(:,8) ~= blockID, :)=[];
        c=1;
        for ii=contrast
          trialsPerContrast(c) = length(data(data(:,2)==ii,2));
          if trialsPerContrast(c) < 0.05*length(data)
          includeContrast(c) = 0;
          end
        c=c+1;
        end

        %% Get psychometric curves conditional to prev/next correct L/R choices
        [ stims, nextChoice, prevChoice ] = Cond_to_Water_PC_Prev_and_Next( dataAll, includeContrast, blockID );
        % Each has 2 rows - after correct and error / Stims only has 1 row

         
        allStimSorted = sort(stims{1},'ascend');
        diffStim = [ allStimSorted(2)  allStimSorted(end-1)];
        for jj=1:length(stims{1})
           if stims{1}(jj) == 0
               stimID = 3;
           elseif stims{1}(jj) > diffStim(2)
               stimID = 5;
           elseif stims{1}(jj) < diffStim(1)
               stimID = 1;
           elseif stims{1}(jj) > 0
               stimID = 4;
           else
               stimID = 2;
           end
%            if any(isnan(shifts(stimID, :, b, AID)))
           if isnan(shifts(stimID, 1, b, AID))
               tmp = nextChoice{1}( :, jj) - prevChoice{1}(:, jj);
%                if ~any(isnan(tmp))
               if any(~isnan(tmp))
               
                   if ~isnan(tmp(1))
                   shifts(stimID, 1, b, AID) = tmp(1);
                   shift_n(stimID,1, b, AID) = 1;
                   end
                   
                   if ~isnan(tmp(2))
                   shifts(stimID, 2, b, AID) = tmp(2);
                   shift_n(stimID,2, b, AID) = 1;
                   end
%                else
%                    shifts(stimID, :, b, AID) = NaN;
               end
           else
               tmp = nextChoice{1}( :, jj) - prevChoice{1}(:, jj);
%                if ~any(isnan(tmp)) 
                if any(~isnan(tmp))  
                   if ~isnan(tmp(1))
                    shifts(stimID, 1, b, AID) = shifts(stimID, 1, b, AID) + tmp(1);
                    shift_n(stimID,1, b, AID) = shift_n(stimID, 1, b, AID)+1;
                   end
                   if ~isnan(tmp(2))
                       if ~isnan( shifts(stimID, 2, b, AID))
                        shifts(stimID, 2, b, AID) = shifts(stimID,2, b, AID) + tmp(2);
                        shift_n(stimID,2, b, AID) = shift_n(stimID,2, b, AID)+1;
                       else
                        shifts(stimID, 2, b, AID) = tmp(2);
                        shift_n(stimID,2, b, AID) = 1;
                       end
                   end
               end
           end
        end
      
        end
    end

    shifts = shifts./shift_n;
    % numRuns = sum(includeAnimal);
    % shifts = shifts(:,:,:, includeAnimal==1);

end