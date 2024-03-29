function [ contrasts, shift] = CondPsychCurves( TStimuli, TBlocks, action, correct, varargin )
%% Returns shifts. Conditional to correct Left/Right actions - with separate arguments
    nTrials = length(TStimuli);
    nIters = size(action,2);
    
    contrasts = unique(TStimuli)';
    blocks = unique( TBlocks)';%[1 2];%
    
    
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
        lst='-';
    end
    % Change responses to 0 and 1
    if any(action == -1 )
        action = (1 + action)/2;
    end
    
    %             Left/Right,   Current stimulus,  Block number
    PC_next = zeros(  2,        length(contrasts), length(blocks), nIters );
    PC_prev = PC_next;
%   shift = nan(  2, length(contrasts), length(blocks) );
       
    b=0;
    for blockID = blocks
    b=b+1; c= 0;
        for stim = contrasts
        c= c+1;
        
        
        for iter = 1:nIters
        % After left correct
            id = [false; TStimuli(2:end) == stim & TBlocks(2:end)==blockID & ...
                          action(1:end-1, iter) == 0 & correct(1:end-1, iter) == 1];
            PC_next(1,c,b,iter) = nanmean( action(id, iter) ); 
            
        % Before left correct
            id = [TStimuli(1:end-1) == stim & TBlocks(1:end-1)==blockID & ...
                   action(2:end, iter) == 0 & correct(2:end, iter) == 1; false];
            PC_prev(1,c,b,iter) = nanmean( action(id, iter) ); 
        
            
        % After right correct
            id = [false; TStimuli(2:end) == stim & TBlocks(2:end)==blockID & ...
                          action(1:end-1, iter) == 1 & correct(1:end-1, iter) == 1];
            PC_next(2,c,b,iter) = nanmean( action(id, iter) ); 
            
        % Before right correct
            id = [TStimuli(1:end-1) == stim & TBlocks(1:end-1)==blockID & ...
                   action(2:end, iter) == 1 & correct(2:end, iter) == 1; false];
            PC_prev(2,c,b,iter) = nanmean( action(id, iter) ); 
        end

        
        end
    end
%     PC_next = PC_next./nIters;
%     PC_prev = PC_prev./nIters;
    
    shift = nanmean(PC_next - PC_prev,4);
    
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
           ic = ~isnan(shift(1,:,b) );
           h(2*b-1) = plot( contrasts(ic), shift(1,ic, b), 'Color', 'b', 'marker', MarkerShape{blockID}, 'markerfacecolor', 'b', 'markersize', 10 , 'linestyle', lst);
           ic = ~isnan(shift(2,:,b) );
           h(2*b  ) = plot( contrasts, shift(2,ic, b), 'Color', 'r', 'marker', MarkerShape{blockID}, 'markerfacecolor', 'r', 'markersize', 10 , 'linestyle', lst);
           
           MLabel{2*b-1} = sprintf('Correct left (Block %s)', BlockLabel{blockID} );
           MLabel{2*b  } = sprintf('Correct right (Block %s)', BlockLabel{blockID} );
       end
       if newfig
         legend(h, MLabel{:}, 'Location', 'SouthEast' ) 
       end
       ylabel('Change in P(R) after correct actions')
       xlabel('Stimulus')
       ylim([ -0.2 0.2])
    end
end