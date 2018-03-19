p = uigetdir('Select file directory of Figures...'); 
cd(p); 
%% Open figures and save as png
F = dir; 

F(1:2) = []; 
delInd = find([F.bytes] == 0); 
F(delInd) = []; 
for i = 1:length(F)
   fig = open(F(i).name); 
   saveas(fig, ['P' F(i).name(1:end-4)], 'png'); 
end

close all