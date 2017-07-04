function [ contrasts, PC] = PsychCurves( TStimuli, TBlocks, action, varargin )
%% Psychometric curves with separate arguments

    nTrials = length(TStimuli);
    nIters = size(action,2);
    
    contrasts = unique(TStimuli)';
    blocks = unique( TBlocks)';
    
    
    if nargin>3
        toplot = varargin{1};
    else
        toplot = false;
    end
    if nargin > 4
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
    
    PC = zeros( length(contrasts), length(blocks) , nIters); 
    IC = zeros(length(contrasts), length(blocks));   
    
    b=0;
    for blockID = blocks
    b=b+1; c= 0;
        
        for stim = contrasts
        c= c+1;
        
        id = ( TStimuli == stim & TBlocks==blockID);
        if ~isempty(id)
            IC( c, b ) = 1;
        end
        for iter = 1:nIters
           PC(c,b, iter) = nanmean( action(id, iter) ); 
        end

        
        end
    end
    PC = nanmean(PC,3);
    
    
    if toplot
       BlockColor = {'b', 'r', 'g', 'k'};
       MarkerShape = {'o', 's', '^', 'v'};
       BlockLabel = {'DA-L', 'DA-R', 'DA-Both', 'DA-None'};
       if newfig
        figure; hold on
       end
       b=0; h = zeros(length(blocks),1);
       for blockID=blocks
           b=b+1;
           h(b) = plot( contrasts(IC(:,b)==1), PC(IC(:,b)==1, b), 'Color', BlockColor{blockID}, 'marker', MarkerShape{blockID}, 'markerfacecolor', BlockColor{blockID}, 'markersize', 10 , 'linestyle', lst);
       end
       if newfig
           legend(h, BlockLabel{blocks}, 'Location', 'SouthEast' ) 
       end
       ylim([ 0 1])
    end
end