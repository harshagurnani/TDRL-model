function PlotColourMap_softmax(alpha,DA_val,noiseSTD,beta,FvalStore)
% function to plot a colour map of goodness of fit determined by
% FvalStoreAlt

% find the location of minimum value of FvalStoreAlt
[~,locMin] = min(FvalStore(:));
[a,b,c,d] = ind2sub(size(FvalStore),locMin);


% plot noiseSTD against beta
colourVal = nan(length(noiseSTD),length(beta));
for j=1:length(beta)
    for i=1:length(noiseSTD)
        
        % Set labels for the parameters we are varying
        xData = 'beta';
        yData = 'noiseSTD';
        
        % Prepare variables for colour map
        x = [min(beta),max(beta)];
        y = [min(noiseSTD),max(noiseSTD)];
        colourVal(i,j) = FvalStore(a,b,i,j);
        
    end
end

% Plot the colour map
figure
imagesc(x,y,colourVal);
    colorbar
    xlabel(xData)
    ylabel(yData)
    dim = [.35 .69 .3 .3];
    annot = 'foo';
    annotation('textbox',dim,'String',annot,'FitBoxToText','on');  


% plot alpha against DA_val
colourVal = nan(length(alpha),length(DA_val));
for j=1:length(DA_val)
    for i=1:length(alpha)
        
        % Set labels for the parameters we are varying
        xData = 'DAval';
        yData = 'alpha';
        
        % Prepare variables for colour map
        x = [min(DA_val),max(DA_val)];
        y = [min(alpha),max(alpha)];
        colourVal(i,j) = FvalStore(i,j,c,d);
        
    end
end

% Plot the colour map
figure
imagesc(x,y,colourVal);
    colorbar
    xlabel(xData)
    ylabel(yData)
    dim = [.35 .69 .3 .3];
    annot = 'foo';
    annotation('textbox',dim,'String',annot,'FitBoxToText','on');

end

