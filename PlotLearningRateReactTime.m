function PlotLearningRateReactTime(Data,data_model)

data_mice=Data.data;

smoothdata_mice(:) = (1 + data_mice(:,3)) / 2;
smoothdata_model(:) = (1 + data_model(:,16)) / 2;
responseTime(:) = data_mice(:,7);

smoothing_factor = 7;
smoothdata_mice(:) = smooth(smoothdata_mice(:),smoothing_factor);
smoothdata_model(:) = smooth(smoothdata_model(:),smoothing_factor);

plotLen = 30;

status = zeros(length(unique(data_mice(:,9))),5);

% we initialse our firstTrials arrays
firstTrials_mice = zeros(plotLen,4);
firstTrials_model = zeros(plotLen,4);
firstTrials_response = zeros(plotLen,4);
% we give these arrays 4 columns to separate out the blocks depending on
% their status, i.e.
% firstTrials(:,1) = 1st session of the day, block 1
% firstTrials(:,2) = 1st session of the day, block 2
% firstTrials(:,3) = 2nd session of the day, block 1
% firstTrials(:,4) = 2nd session of the day, block 2

c=1;
counter = zeros(1,4);
for blockN = unique(data_mice(:,9))'
   
   day_session = data_mice(data_mice(:,9)==blockN,11);
   blockID = data_mice(data_mice(:,9)==blockN,8);
   
   % 'status' columns are as follows
   % col 1: block number as given in the data
   % col 2: number of trials in the block
   % col 3: 1 for first session of the day, 2 otherwise
   % col 4: block ID (which defines reward structure) (1 or 2)
   % col 5: block ID of previous block (0 for first block)
   status(c,:) = [blockN, length(blockID), min(day_session(1),2),...
      blockID(1), status(max(1,c-1),4)];
   
   
   if status(c,2)<plotLen
      c=c+1;
      continue
   end
   
   if status(c,4)==status(c,5) && status(c,3)~=1
      % if the blockID and the previous block ID are the same, and it isn't
      % the first session of the day, skip this block
      c=c+1;
      continue
   end
   
   
   tmp_mice = smoothdata_mice(data_mice(:,9)==blockN);
   tmp_model = smoothdata_model(data_mice(:,9)==blockN);
   tmp_response = responseTime(data_mice(:,9)==blockN);
   
   firstTrials_mice(:,(status(c,3)-1)*2 + status(c,4)) = ...
      firstTrials_mice(:,(status(c,3)-1)*2 + status(c,4)) ...
      + tmp_mice(1:plotLen)';
   firstTrials_model(:,(status(c,3)-1)*2 + status(c,4)) = ...
      firstTrials_model(:,(status(c,3)-1)*2 + status(c,4)) ...
      + tmp_model(1:plotLen)';
   firstTrials_response(:,(status(c,3)-1)*2 + status(c,4)) = ...
      firstTrials_response(:,(status(c,3)-1)*2 + status(c,4)) ...
      + tmp_response(1:plotLen)';
   
   counter((status(c,3)-1)*2 + status(c,4)) = ...
      counter((status(c,3)-1)*2 + status(c,4)) + 1;

   c=c+1;

end

plt_mice = firstTrials_mice ./ counter;
plt_model = firstTrials_model ./ counter;


plt_response(:,1) = ...
   (firstTrials_response(:,1) + firstTrials_response(:,2)) / ...
   (counter(1) + counter(2));
plt_response(:,2) = ...
   (firstTrials_response(:,3) + firstTrials_response(:,4)) / ...
   (counter(3) + counter(4));



headingAction = {sprintf(['Average behaviour early in a block:\n',...
   'First block of the day']),...
   sprintf(['Average behaviour early in a block:\n',...
   'Not first block of the day'])};

headingResponse = {sprintf(['Average reaction time early in a block:\n',...
   'First block of the day']),...
   sprintf(['Average reaction time early in a block:\n',...
   'Not first block of the day'])};

fig = figure;

actionAx = zeros(2,1);
actionAx(1) = axes(fig);
actionAx(2) = axes(fig);
responseAx = zeros(2,1);
responseAx(1) = axes(fig);
responseAx(2) = axes(fig);

for i=1:2

   actionAx(i)=subplot(2,2,i); hold on;
   plot(plt_mice(:,(i-1)*2 + 1),'color','b','linewidth',1.2)
   plot(plt_model(:,(i-1)*2 + 1), 'color','b','linestyle','--',...
      'linewidth',1.2)
   plot(plt_mice(:,(i-1)*2 + 2),'color','r','linewidth',1.2)
   plot(plt_model(:,(i-1)*2 + 2), 'color','r','linestyle','--',...
      'linewidth',1.2)
   ylim([0 1])
   
   ylabel('Fraction rightward choice')
   xlabel('Number of trials after block start')
   
   title(actionAx(i),headingAction(i))
   
   responseAx(i)=subplot(2,2,2+i); hold on;
   plot(plt_response(:,i),'color','k','linewidth',1.2)
   
   ylim([0 max(plt_response(:))+1])
   ylabel(sprintf('Normalised average\nresponse time (s)'))
   xlabel('Number of trials after block start')
   
   title(responseAx(i),headingResponse(i))
      
end

legend(actionAx(1),'Blk 1 mouse','Blk 1 model','Blk 2 mouse',...
   'Blk 2 model','Location','northwest')

%legend(responseAx(2),'Blk 1 response','Blk 2 response','Location','north')


end