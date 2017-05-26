 function AllThePlots(Data,data_model, xanswerMin, varargin)
% Combining all the plots we have

data_mice = Data.data;
animalID = Data.ID;

% changing y-axis values
data_mice(:,3) = (1 + data_mice(:,3)) ./ 2;
data_model(:,16) = (1 + data_model(:,16)) ./ 2;
data_model(:,22) = (1 + data_model(:,22)) ./ 2;

contrasts = unique(data_mice(:,2))';
blocks = unique(data_mice(:,8))';

if nargin>3
   action = varargin{1};
else
   action = data_model(:,22);
end

if nargin>4
   correct = varargin{2};
else
   correct = 1*(action==data_mice(:,3) & data_mice(:,10)==1) + 1*(action~=data_mice(:,3) & data_mice(:,10)==0);
end

if nargin>5
    includeContrast=varargin{3};
else
    includeContrast = ones(length(contrasts)*length(blocks),1);
end
ic = reshape(includeContrast, length(contrasts), length(blocks)
color=[0 0 1; 1 0 0; 0 1 0; 0 0 0];   % blue, red, green, black

full = figure;

%% PLOT 1 - Psychometric curves
f1=subplot(3,3,1);
PlotPC_Mice_Model(data_mice, data_model, action, includeContrast);%, full, [3 3 1] )

%% PLOT 2 
%% Conditional curve after Left
f2=subplot(3,3,2);
PlotCondPC_Mice_Model(data_mice, data_model, action, correct,  includeContrast, [1 2]);%, full, [3 3 2] )

%% PLOT 3 
%% Conditional curve after right
f3=subplot(3,3,3);
PlotCondPC_Mice_Model(data_mice, data_model, action, correct,  includeContrast, [3 4]);%, full, [3 3 3] )

% %% PLOT 2
% 
% plotLen = length(data_mice);
% if length(data_mice)<plotLen
%    
%    plotLen = length(data_mice);
%    
% end
% plotLenVec = cat(1,ones(plotLen,1),zeros(length(data_mice)-plotLen,1));
% 
% mouse = data_mice(plotLenVec==1,3) ;
% model = data_model(plotLenVec==1,22);
% 
% smoothmouse = mouse*0.96 + 0.02;
% smoothmodel = model*0.96 + 0.02;
% 
% smoothmouse = smooth(smoothmouse,10);
% smoothmodel = smooth(smoothmodel,10);
% 
% 
% blockID = data_mice(plotLenVec==1,8) - 1;
% blockID(blockID==1) = 1.03;
% blockID(blockID==0) = -0.03;
% 
% 
% subplot(3,1,2,'Parent',full);  hold on
% plot(mouse,'marker','x','MarkerEdgeColor',[0.2 0 0.6],'markerfacecolor',...
%    [0.2 0 0.6],'markersize',2,'linestyle','none')
% 
% plot(blockID,'marker','o','MarkerEdgeColor',[1 .5 0],'markerfacecolor',...
%    [1 .5 0],'markersize',2,'linestyle','none')
% 
% plot(smoothmouse,'color',[0.2 0 0.6])
% 
% %plot(smoothmodel,'color',[0.7 0.7 0.7])
% plot(smoothmodel,'color',[0 0 0])
% 
% 
% ylim([-0.05 1.05])
% xlim([0 plotLen])
% 
% xlabel('Trial number')
% ylabel('Fraction rightward choice')
% 
% %pbaspect([4 1 1])
% set(gca,'visible','on');
% 
% 
% %% PLOT 3
% 
% smoothing_factor = 11;
% smoothdata_mice(:) = smooth(data_mice(:,3),smoothing_factor);
% smoothdata_model(:) = smooth(data_model(:,16),smoothing_factor);
% 
% plotLen = 30;
% 
% status = zeros(length(unique(data_mice(:,9))),5);
% 
% firstTrials_mice = zeros(plotLen,4);
% firstTrials_model = zeros(plotLen,4);
% 
% c=1;
% counter = zeros(1,4);
% for blockN = unique(data_mice(:,9))'
%    
%    day_session = data_mice(data_mice(:,9)==blockN,11);
%    blockID = data_mice(data_mice(:,9)==blockN,8);
%    
%    % 'status' columns are as follows
%    % col 1: block number as given in the data
%    % col 2: number of trials in the block
%    % col 3: 1 for first session of the day, 2 otherwise
%    % col 4: block ID (which defines reward structure) (1 or 2)
%    % col 5: block ID of previous block (0 for first block)
%    status(c,:) = [blockN, length(blockID), min(day_session(1),2),...
%       blockID(1), status(max(1,c-1),4)];
%    
%    
%    if status(c,2)<plotLen
%       c=c+1;
%       continue
%    end
%    
%    if status(c,4)==status(c,5) && status(c,3)~=1
%       % if the blockID and the previous block ID are the same, and it isn't
%       % the first session of the day, skip this block
%       c=c+1;
%       continue
%    end
%    
%    
%    tmp_mice = smoothdata_mice(data_mice(:,9)==blockN);
%    tmp_model = smoothdata_model(data_mice(:,9)==blockN);
%    
%    firstTrials_mice(:,(status(c,3)-1)*2 + status(c,4)) = ...
%       firstTrials_mice(:,(status(c,3)-1)*2 + status(c,4)) ...
%       + tmp_mice(1:plotLen)';
%    firstTrials_model(:,(status(c,3)-1)*2 + status(c,4)) = ...
%       firstTrials_model(:,(status(c,3)-1)*2 + status(c,4)) ...
%       + tmp_model(1:plotLen)';
%    
%    counter((status(c,3)-1)*2 + status(c,4)) = ...
%       counter((status(c,3)-1)*2 + status(c,4)) + 1;
% 
%    c=c+1;
% 
% end
% 
% plt_mice = firstTrials_mice ./ counter;
% plt_model = firstTrials_model ./ counter;
% 
% 
% heading = {sprintf(['Average behaviour early in a block:\n',...
%    'First block of the day']),...
%    sprintf(['Average behaviour early in a block:\n',...
%    'Not first block of the day'])};
% 
% 
% for i=1:2
% 
%    subplot(3,2,4+i,'Parent',full); hold on
%    plot(plt_mice(:,(i-1)*2 + 1),'color','b','linewidth',1.2)
%    plot(plt_model(:,(i-1)*2 + 1), 'color','b','linestyle','--',...
%       'linewidth',1.2)
%    plot(plt_mice(:,(i-1)*2 + 2),'color','r','linewidth',1.2)
%    plot(plt_model(:,(i-1)*2 + 2), 'color','r','linestyle','--',...
%       'linewidth',1.2)
%    ylim([0 1])
%    
%    ylabel('Fraction rightward choice')
%    xlabel('Number of trials after block start')
%    
%    title(heading(i))
% 
% end
%    
%    legend('Blk 1 mouse','Blk 1 model','Blk 2 mouse','Blk 2 model',...
%    'Location','north')

end

