function [ L_T ] = RemoveT_LR( T )
%REMOVET_LR Summary of this function goes here
%   Function that removes all but the first two blocks (to get rid of learning effect)

subs = unique(T{:,1});

for sub = 1 : length(subs)
    
    for stim = 1 : 2
        
        sset = T{:,1} == sub & T{:,2} == stim & T{:,3} == 1; 
        maxBlocks = max(T{sset,4}); 
        T(sset & ~ismember(T{:,4},[maxBlocks,maxBlocks-1]),:) = []; 
        
    end
    
end

L_T = T; 

end

