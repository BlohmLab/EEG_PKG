function [ mu_rxn_samp ] = GetMedianRXN( E )
%GETMEDIANRXN Get median reaction time of EEG structure

disp('Test start'); 

X = [E(:).Anodal E(:).Cathodal]; 

X1 = squeeze(struct2cell(X)); 
Y = X1(:);
Z = [Y{:}];
T = cat(3,Z(:).T); 
RXN = squeeze(T(:,3,:)); 
baseRXN = squeeze(T(:,1,:)); 
RXN = RXN-baseRXN;
RXN = nanmedian(RXN,1);
mu_rxn_samp = nanmean(RXN); 

end

