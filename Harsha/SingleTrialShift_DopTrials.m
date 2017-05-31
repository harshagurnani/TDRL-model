function [ shifts ] = SingleTrialShift_DopTrials( dataGiven , varargin )

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
    includeContrast = ones(length(contrast),1);

    
    
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
        
        b=0;
        for blockID = blocks
        b=b+1;
        %% Get contrats to include
        data=dataAll; 
        data( data(:,2) ~= blockID, :)=[];
        c=1;
        for ii=contrast
          trialsPerContrast(c) = length(data(data(:,2)==ii,2));
          if trialsPerContrast(c) < 0.05*length(data)
          includeContrast(c) = 0;
          end
        c=c+1;
        end

        %% Get psychometrical curves conditional to prev/next correct L/R choices
        [ stims, nextChoice, prevChoice ] = Cond_to_Dop_PC_Prev_and_Next( dataAll, includeContrast, blockID );
        % Each has 2 rows - after correct and error / Stims only has 1 row

        diffStim = [ max(stims{1}(stims{1}<0))  min(stims{1}(stims{1}>0))];
        for jj=1:length(stims{1})
           if stims{1}(jj) == 0
               stimID = 3;
           elseif stims{1}(jj) == diffStim(1)
               stimID = 2;
           elseif stims{1}(jj) == diffStim(2)
               stimID = 4;
           elseif stims{1}(jj) > diffStim(2)
               stimID = 5;
           else
               stimID = 1;
           end
           if any(isnan(shifts(stimID, :, b, AID)))
               tmp = nextChoice{1}( :, jj) - prevChoice{1}(:, jj);
               if ~any(isnan(tmp))
                   shifts(stimID, :, b, AID) = tmp';
                   shift_n(stimID,:, b, AID) = 1;
               else
                   shifts(stimID, :, b, AID) = NaN;
               end
           else
               tmp = nextChoice{1}( :, jj) - prevChoice{1}(:, jj);
               if ~any(isnan(tmp))  
                   shifts(stimID, :, b, AID) = shifts(stimID, :, b, AID) + tmp';
                   shift_n(stimID,:, b, AID) = shift_n(stimID, :, b, AID)+1;
               end
           end
        end
      
        end
    end

    shifts = shifts./shift_n;
    % numRuns = sum(includeAnimal);
    % shifts = shifts(:,:,:, includeAnimal==1);

end