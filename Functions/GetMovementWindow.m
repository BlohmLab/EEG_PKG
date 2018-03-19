function [ nineP ] = GetMovementWindow( E )
%GETMOVEMENTWINDOW Summary of this function goes here
%   Detailed explanation goes here
X = [E(:).Anodal E(:).Cathodal]; 

X1 = squeeze(struct2cell(X)); 
Y = X1(:);
Z = [Y{:}];
T = cat(3,Z(:).T); 
RXN = squeeze(T(:,3,:)); 
baseRXN = squeeze(T(:,2,:)); 
RXN = RXN-baseRXN;
RXN = RXN(:); 
RXN(isnan(RXN)) = []; 
mRXN = median(RXN(:)); 

%Get 90th percentile of movements
q = [0.05 0.95]; 
nineP = quantile(RXN,q); 



end

