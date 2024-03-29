function [ stimuli, shifts ] = CondPsychCurves_StimLevel_and_WaterCorrect( TStimuli, TBlocks, action, correct, varargin )
%% Shift in P(R) after water reward on easy versus difficult stimuli

    nTrials = length(TStimuli);
    nIters = size(action,2);

    contrasts = unique(TStimuli)';
    blocks = [1 2];


    if nargin>4
        toplot = varargin{1};
    else
        toplot = false;
    end
    if nargin > 5
        newfig = varargin{2};
        lst = '--';
    else
        newfig = true;
        lst='--';
    end
    
    % Change responses to 0 and 1
    if any(action == -1 )
        action = (1 + action)/2;
    end

    % Easy/difficult stimulus, current stimulus, blockID
    PC_next = zeros( 2, length(contrasts), length(blocks), nIters );
    PC_prev = PC_next;
    
    WaterChoice = {1, 0};
    WaterGiven = nan(length(nTrials),1);


    markerShape = {'o', 's', '^', 'v'};
    blockColor = {'b','r','g','k'};
    blockLabel = {'Left', 'Right', 'DA-Both', 'DA-None'};


    includeContrast = ones(length(contrasts),1);
    trialsPerContrast = nan(length(contrasts), 1);
    c=1;
    for ii=contrasts
        trialsPerContrast(c) = length(TStimuli==ii);
        if trialsPerContrast(c) < 0.05*nTrials
          includeContrast(c) = 0;
        end
        c=c+1;
    end

    stimuli = contrasts(includeContrast ==1 );
    ic = includeContrast==1;
    shifts   = nan(  2, length(stimuli), length(blocks), nIters );
    

    IsStimEasy = nan(nTrials,1);
    IsStimDiff = IsStimEasy;

    easystim = [ min(stimuli)              max(stimuli) ];
    diffstim = [ max(stimuli(stimuli<0)) 0 min(stimuli(stimuli>0)) ];

    IsStimEasy = (TStimuli<=easystim(1) | TStimuli>=easystim(2));
    IsStimDiff = (TStimuli>=diffstim(1) & TStimuli<=diffstim(2));

    IsStimEasy = reshape(IsStimEasy, [nTrials,1]);
    IsStimDiff = reshape(IsStimDiff, [nTrials,1]);

    for iter = 1:nIters
    b=0;
    for blockID = blocks
       b=b+1;

       WaterGiven = (action(:,iter)==WaterChoice{blockID});
       WaterGiven = reshape(WaterGiven, [nTrials, 1]);

       % Indx 1 is after easy,indx 2 is after diff stim.
       c=1; %current stimuli
       for istim = contrasts

          %After/Before easy stim+dopamine.

          id=[ false; TStimuli(2:end)==istim  & TBlocks(2:end)==blockID & ...           % current stimulus
              WaterGiven(1:end-1,1) & correct(1:end-1,iter)==1 & IsStimEasy(1:end-1,1)]; % previous action
          PC_next(1,c,blockID,iter)= nanmean(action(id,iter)) ;
          id=[ TStimuli(1:end-1)==istim  & TBlocks(1:end-1)==blockID & ...  
              WaterGiven(2:end,1) & correct(2:end,iter)==1     & IsStimEasy(2:end); false];
          PC_prev(1,c,blockID,iter)= nanmean(action(id,iter)) ;

          %After/Before difficult stim + dopamine
          id=[false; TStimuli(2:end)==istim  & TBlocks(2:end)==blockID & ...           % current stimulus
              WaterGiven(1:end-1,1) & correct(1:end-1,iter)==1 & IsStimDiff(1:end-1)]; % previous action
          PC_next(2,c,blockID,iter)= nanmean(action(id,iter)) ;
          id=[ TStimuli(1:end-1)==istim  & TBlocks(1:end-1)==blockID & ...  
              WaterGiven(2:end,1) & correct(2:end,iter)==1     & IsStimDiff(2:end); false];
          PC_prev(2,c,blockID,iter)= nanmean(action(id,iter)) ;


          c=c+1;
       end


       

   end
   end
   shifts = nanmean(PC_next(:,ic,:,:) - PC_prev(:,ic,:,:),4);
   
   if toplot
        if newfig
        figure; 
        end
        hold on;

        h = zeros(4,1);
        mLabel = cell(4,1);

        b=0;
        for blockID=blocks
        b=b+1;
        h(2*b-1) = plot( stimuli, nanmean(PC_next(1,ic,blockID,:),4), 'marker', markerShape{blockID},...
                    'color', blockColor{blockID}, 'linestyle','--');

        h(2*b) = plot( stimuli, nanmean(PC_next(2,ic,blockID,:),4), 'marker', markerShape{blockID},...
                    'color', blockColor{blockID}, 'linestyle','-','markerfacecolor', blockColor{blockID});
        mLabel{2*b-1} = sprintf(' %s - After Easy Stim', blockLabel{blockID});                          
        mLabel{2*b}   = sprintf(' %s - After Diff Stim', blockLabel{blockID});                          
        end
        legend(h, mLabel{:})
        ylabel('P(R)')
        xlabel('Stimulus')
   end

end

