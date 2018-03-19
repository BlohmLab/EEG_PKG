function [ linData ] = LinearizePhase( circData )
%LINEARIZEPHASE Take circular data and linearizes it by appending jumps due
%to circular domain representation 

%Determine jump thresholds 
jmpThres = pi; 

%Determine possible sampling error from pi 
%Update this to reflect phase trajectory velocity and sampling error
%interpolation
epsilon = 0.2; 

%Try looped frequency version
[freqSiz,~] = size(circData); 
linData = circData; 

%For each frequency
for freq = 1 : freqSiz 
    
    %Get current frequency
    freqDat = circData(freq,:); 
    
    %Compute difference
    dtF = diff(freqDat); 
    
    %Check if discontinuity exists and determine directionality 
    jmpInd = find(abs(dtF) > (jmpThres - epsilon)); 
    
    %For each jump append rest of time-series 
    for jmp = 1 : length(jmpInd)
        
        linData(freq,jmpInd(jmp)+1:end) = linData(freq,jmpInd(jmp)+1:end) - dtF(jmpInd(jmp)); 
        
    end
    
    
    
end


end

