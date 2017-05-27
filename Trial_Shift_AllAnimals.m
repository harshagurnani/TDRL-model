blocks = [1 2];

allAnimals = dir('*.mat');
numAnimals = length(allAnimals);

stimuli = [-1 -0.75 -0.5 -0.4 -0.25 -0.1  -0.05 0 0.05  0.1 0.25 0.4 0.5 0.75 1];
shifts = nan( length(stimuli), 2, length(blocks), numAnimals );
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
    data( data(:,2) ~= blockID, :)=[];
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
    [ stims, nextChoice, prevChoice ] = CondPC_Prev_and_Next( dataAll, includeContrast, blockID );
    % Each has 2 rows - after left and right
    for jj=1:length(stims{1})
       minDist = min( abs(stims{1}(jj)-stimuli));
       stimID = find( abs(stims{1}(jj)-stimuli) == minDist );
       shifts(stimID, :, b, AID) = nextChoice{1}( :, jj) - prevChoice{1}(:, jj);
    end
    
    end
end