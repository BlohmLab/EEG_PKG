function [ bRegMat ] = SortBehavReg( data, grp )
%BEHAVREGSORT Summary of this function goes here
%   Take in a grouping value and return a behavioural matrix for performing
%   regressions
%For subtractive comparisons standard is L - R (match with EEG)

c_group = {'TIHP','targIHP','reachVec','EHemifield',...
    'LHemifield','RHemifield','Reach','Target','IHP'};

subNum = unique(data(:,5)); 

%Run conditionals for grouping types and check
if strcmp(grp,c_group{2}) %Special case (perform regression across targets for each initial hand position)
    
    %Compute a slope for each subject for the regression
    for s = 1 : length(subNum)
        
        %Select subject
        subData = data(data(:,5) == s,:); 
        
        %Select Left IHP
        L = subData(subData(:,3) == -7.5,:);
        R = subData(subData(:,3) == 7.5,:);
        
        %Perform regression for this subject 
        regL = regress(L(:,1), [ones(size(L,1),1) L(:,2)]); 
        regR = regress(R(:,2), [ones(size(R,1),1) R(:,2)]);
        
        %Store subject difference regression 
        bRegMat(s) = regL(2) - regR(2); 
        
    end
    
elseif strcmp(grp,c_group{3}) %Special case (perform reach type regression) 
    
    %Treat duplicates as independent values (don't average them) 
    
    for s = 1 : length(subNum) 
        
        subData = data(data(:,5) == s,:); 
        reg = regress(subData(:,1), [ones(size(subData(:,1),1),1) subData(:,8)]); 
        bRegMat(s) = reg(2); 
        
    end
    
elseif strcmp(grp,c_group{4:end}) %General subtraction case
    
    for s = 1 : length(subNum)
        
        subData = data(data(:,5) == s,:);
        
        %Now specific cases, get subtraction matrices (update this to
        %reflect entire case grouping for horizontal error computation) 
        if strcmp(grp,'EHemifield')
            
            L = subData(subData(:,2) < 0 & subData(:,3) < 0,:); 
            R = subData(subData(:,2) > 0 & subData(:,3) > 0,:); 
            
        elseif strcmp(grp,'RHemifield')
            
            L = subData(subData(:,3) > 0 & subData(:,2) > 0,:); 
            R = subData(subData(:,3) > 0 & subData(:,2) < 0,:); 
            
        elseif strcmp(grp,'LHemifield')
            
            L = subData(subData(:,3) < 0 & subData(:,2) < 0,:); 
            R = subData(subData(:,3) < 0 & subData(:,2) > 0,:); 
            
        elseif strcmp(grp,'Reach')
            
            L = subData(subData(:,8) < 0,:); 
            R = subData(subData(:,8) > 0,:); 
            
        elseif strcmp(grp,'Target')
            
            L = subData(subData(:,2) < 0,:); 
            R = subData(subData(:,2) > 0,:); 
            
        elseif strcmp(grp,'IHP')
            
            L = subData(subData(:,3) < 0,:); 
            R = subData(subData(:,3) > 0,:); 
            
        end
        
        bRegMat(s) = L - R; 
        
    end
    
end



end

