blocks = [1 2];

allAnimals = dir('*.mat');
numAnimals = length(allAnimals);

% 2 = Easy , 1 = Difficult (second lowest stimulus), 0 = zero contrast
stimuli = [ -2  -1  0  1  2 ];
shifts = nan( length(stimuli), 2, length(blocks), numAnimals );
shift_n = nan(length(stimuli),2, length(blocks), numAnimals);
includeAnimal = zeros(numAnimals,1);
% stimulus, conditional to left/right actions, block, animalID
for AID=1:numAnimals
    dataAll = importdata(allAnimals(AID).name);
    contrast = unique(dataAll(:,2))';
    b=0;
    for blockID = blocks
    b=b+1;
    %% Get contrats to include
    % All contrasts in data
    data=dataAll; 
    data( data(:,8) ~= blockID, :)=[];
    trialsPerContrast = nan(length(contrast),1);
    includeContrast = ones(length(contrast),1);
   
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
    if ~any(abs(stims{1}) > 0.6)
        includeAnimal(AID) = 1;
        % Not include animals with contrast levels > 0.6
        stimWithoutZero = stims{1}(stims{1}~=0);

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
end

shifts = shifts./shift_n;
shifts_backup = shifts;

% numAnimals = sum(includeAnimal);
% shifts = shifts(:,:,:, includeAnimal==1);