function [ polDiff ] = computeDiffTFR( dSTFR )
%COMPUTEDIFFTFR
%   House-keeping function to ensure that TFR analysis is suited for plotting
%   If left/right analysis is wanted will take difference between left and
%   rightward saccade motion
%   If not, will attempt to squeeze dSTFR to right dimensions
LEFT = 1;
RIGHT = 2;

%Check for dimensions
siz = size(dSTFR,3); 

if (siz == 1)
    polDiff = squeeze(dSTFR); 
else
    polDiff = dSTFR(:,:,LEFT,:) - dSTFR(:,:,RIGHT,:); 
end

end

