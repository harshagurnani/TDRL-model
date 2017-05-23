function PlotColourMap_NLL(alpha,DA_val,noiseSTD,FvalStore, varargin)
% function to plot a colour map of goodness of fit determined by
% FvalStoreAlt

% find the location of minimum value of FvalStore
[~,locMin] = min(FvalStore(:));
[a,b,c] = ind2sub(size(FvalStore),locMin);


% plot alpha against DA_val
clear colourVal x y
colourVal = nan(length(alpha),length(DA_val));
for j=1:length(DA_val)
   for i=1:length(alpha)
      
      % Set labels for the parameters we are varying
      xData = 'DAval';
      yData = 'alpha';
      
      % Prepare variables for colour map
      x = [min(DA_val),max(DA_val)];
      y = [min(alpha),max(alpha)];
      colourVal(i,j) = FvalStore(i,j,c);
      
   end
end

% Plot the colour map
figure('position',[100,100,1200,400]);
if nargin>4
    suptitle(varargin{1})
end
subplot(1,3,1);
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)


% plot alpha against noiseSTD
clear colourVal x y
colourVal = nan(length(alpha),length(noiseSTD));
for j=1:length(noiseSTD)
   for i=1:length(alpha)
      
      % Set labels for the parameters we are varying
      xData = 'noiseSTD';
      yData = 'alpha';
      
      % Prepare variables for colour map
      x = [min(noiseSTD),max(noiseSTD)];
      y = [min(alpha),max(alpha)];
      colourVal(i,j) = FvalStore(i,b,j);
      
   end
end

% plot the colour map
subplot(1,3,2)
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)


% plot DA_val against noiseSTD
clear colourVal x y
colourVal = nan(length(DA_val),length(noiseSTD));
for j=1:length(noiseSTD)
   for i=1:length(DA_val)
      
      % Set labels for the parameters we are varying
      xData = 'noiseSTD';
      yData = 'DAval';
      
      % Prepare variables for colour map
      x = [min(noiseSTD),max(noiseSTD)];
      y = [min(DA_val),max(DA_val)];
      colourVal(i,j) = FvalStore(a,i,j);
      
   end
end

% plot the colour map
subplot(1,3,3)
imagesc(x,y,colourVal);
colorbar
xlabel(xData)
ylabel(yData)


end

