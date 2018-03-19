function PlotBehavCoef( P, varargin )
%PLOTBEHAVCOEF Summary of this function goes here
%   Run correlation, plot then save figure if option enabled

 l = length(P.EEGMat); 
%Run correlation and get R and p value for separate components 
[rA, pValA] = corr(P.EEGMat(1:l/2),P.BMat(1:l/2));
[rC, pValC] = corr(P.EEGMat((l/2)+1:l),P.BMat((l/2)+1:l)); 
[r, pVal] = corr(P.EEGMat,P.BMat); 

%Fit the data
coeffsA = polyfit(P.EEGMat(1:l/2),P.BMat(1:l/2),1);
coeffsC = polyfit(P.EEGMat((l/2)+1:l),P.BMat((l/2)+1:l),1); 
coeffs = polyfit(P.EEGMat,P.BMat,1); 

%Get fitted values
fittedXA = linspace(min(P.EEGMat(1:l/2)), max(P.EEGMat(1:l/2)), 200);
fittedYA = polyval(coeffsA,fittedXA);

fittedXC = linspace(min(P.EEGMat((l/2)+1:l)), max(P.EEGMat((l/2)+1:l)), 200); 
fittedYC = polyval(coeffsC,fittedXC); 

fittedX = linspace(min(P.EEGMat),max(P.EEGMat),200); 
fittedY = polyval(coeffs,fittedX); 

%% Parsing save option and figure show option
p = inputParser;

defaultSave = '0';
checkSave = @(x) isa(x,'char');

defaultPlot = false;
checkPlot = @(x) islogical(x);

addOptional(p,'Save',defaultSave,checkSave);
addOptional(p,'Plot',defaultPlot,checkPlot);

p.KeepUnmatched = true;
parse(p,varargin{:});

%% Plot Option

if (p.Results.Plot)
    
   
    %Plot the data
    figure;
    %Separate out data into two components 
    plot(P.EEGMat(1:l/2),P.BMat(1:l/2),'k.', 'MarkerSize',10); hold on;%Anodal
    plot(P.EEGMat((l/2 + 1):end),P.BMat((l/2 + 1):end),'r.','MarkerSize',10); %Cathodal
    
%     %Run subject sub-plot
%     for i = 1 : l/2
%         plot([P.EEGMat(i) P.EEGMat(i + l/2)], [P.BMat(i) P.BMat(i + l/2)],'--'); hold on;
%     end
    
    %Plot fitted line
    plot(fittedXA,fittedYA,'k--','LineWidth',1); hold on; 
    plot(fittedXC,fittedYC,'r--','LineWidth',1); hold on; 
    plot(fittedX,fittedY,'b-','LineWidth',1); hold on; 
    
    %Set axes and title
    xlabel(['EEG Amplitude Window: ' P.window]); 
    ylabel([P.error ' Error']); 
    title([P.area ' Stim: ' P.stim ' Cond: ' P.cond  ' Group: ' P.group]); 
    
    %Add R^2 and p value to plot
    str = {['r = ' num2str(r)],['R^2 = ' num2str(r^2)], ['p = ' num2str(pVal)]};
    text(0.75,0.75,str,'Units','Normalized');
    
    str = {['Ar = ' num2str(rA)],['AR^2 = ' num2str(rA^2)], ['p = ' num2str(pValA)]};
    text(0.15,0.15,str,'Units','Normalized');
    
    str = {['Cr = ' num2str(rC)],['CR^2 = ' num2str(rC^2)], ['p = ' num2str(pValC)]};
    text(0.15,0.75,str,'Units','Normalized');
       
    legend('Anodal','Cathodal'); 
end

%% Save Option (if save then check if directory added) 
if (~strcmpi(p.Results.Save,'0'))   
    
    %If significant then add an additional tag
    if (pVal < 0.05) 
        tag = '_S'; 
    else
        tag = ''; 
    end
    
    
    saveas(gcf,[p.Results.Save tag],'fig'); 
    saveas(gcf,[p.Results.Save tag],'png'); 
    
end


end

