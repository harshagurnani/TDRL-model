function [ contrasts, shift] = CondPsychCurves_AfterWater( TStimuli, TBlocks, action, correct, varargin )

    nTrials = length(TStimuli);
    nIters = size(action,2);
    
    contrasts = unique(TStimuli)';
    blocks = unique( TBlocks)';
    
    
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
    
    PC_next = zeros(  2, length(contrasts), length(blocks) );
    PC_prev = PC_next;
    shift = nan(  2, length(contrasts), length(blocks) );
    DopChoice = {0, 1, [0 1 ], 2};
    DOPGiven = nan(length(nTrials),1);
    b=0;
    for blockID = blocks
        b=b+1; c= 0;
        DChoice = DopChoice{blockID}
        for stim = contrasts
        c= c+1;
        
        
        for iter = 1:nIters
            for jj=1:nTrials
                DOPGiven(jj,1) = any(action(jj,iter)==DChoice);
            end
        % After Dop Correct
            id = [false; TStimuli(2:end) == stim & TBlocks(2:end)==blockID & ...
                          DOPGiven(1:end-1,1) & correct(1:end-1, iter) == 1];
            PC_next(1,c,b) = PC_next(1,c,b) + nanmean( action(id, iter) ); 
            
        % Before Dop correct
            id = [TStimuli(1:end-1) == stim & TBlocks(1:end-1)==blockID & ...
                   DOPGiven(2:end,1) & correct(2:end, iter) == 1; false];
            PC_prev(1,c,b) = PC_prev(1,c,b) + nanmean( action(id, iter) ); 
        
            
        % After Dop Error
            id = [false; TStimuli(2:end) == stim & TBlocks(2:end)==blockID & ...
                          DOPGiven(1:end-1,1) & correct(1:end-1, iter) == 0];
            PC_next(2,c,b) = PC_next(2,c,b) + nanmean( action(id, iter) ); 
            
        % Before Dop Error
            id = [TStimuli(1:end-1) == stim & TBlocks(1:end-1)==blockID & ...
                   DOPGiven(2:end,1) & correct(2:end, iter) == 0; false];
            PC_prev(2,c,b) = PC_prev(2,c,b) + nanmean( action(id, iter) ); 
        end

        
        end
    end
    PC_next = PC_next./nIters;
    PC_prev = PC_prev./nIters;
    
    shift = PC_next - PC_prev;
    
    if toplot
%        BlockColor = {'b', 'r', 'g', 'k'};
%        lst = '-';
       MarkerShape = {'o', 's', '^', 'v'};
       BlockLabel = {'DA-L', 'DA-R', 'DA-Both', 'DA-None'};
%      
       if newfig
        figure; hold on
       end
       b=0; h = zeros(2*length(blocks),1);
       MLabel = cell(2*length(blocks), 1);
       for blockID=blocks
           b=b+1;
           h(2*b-1) = plot( contrasts, shift(1,:, b), 'Color', [0 0.5 0], 'marker', MarkerShape{blockID}, 'markerfacecolor', [0 0.5 0], 'markersize', 10 , 'linestyle', lst);
           h(2*b  ) = plot( contrasts, shift(2,:, b), 'Color', 'r', 'marker', MarkerShape{blockID}, 'markerfacecolor', 'r', 'markersize', 10 , 'linestyle', lst);
           
           MLabel{2*b-1} = sprintf('After Water correct (Block %s)', BlockLabel{blockID} );
           MLabel{2*b  } = sprintf('After Water error (Block %s)', BlockLabel{blockID} );
       end
       if newfig
         legend(h, MLabel{:}, 'Location', 'SouthEast' ) 
       end
       ylabel('Change in P(R)')
       xlabel('Stimulus')
       ylim([ -0.2 0.2])
    end
end