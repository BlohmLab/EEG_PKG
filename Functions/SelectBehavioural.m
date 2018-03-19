function [ mu, slope ] = SelectBehavioural( varStr, type, condition, stimulation, area )
%SELECTBEHAVIOURAL Summary of this function goes here
%   Detailed explanation goes here

%Column constants
POL = 4;
COND = 7;
SUB = 5;
TARGET = 2;
IHAND = 3;
SITE = 6;
REACH = 8;
%Correlation Constants
MAX_PTHRES = 0.05;
BETA_RANGE = 5;


c_var = {'abs' 'amp' 'x' 'y' 'ddm' 'var'};
c_type = {'Target', 'Reach','LHemifield','IHP', 'EHemifield','TIHP'};

varFind = strfind(c_var,varStr);
varInd = ~cellfun('isempty',varFind);
typeFind = strfind(c_type,type);
typeInd = ~cellfun('isempty',typeFind);
gType = find(typeInd);

X = load('7colvs_v6.mat');

mVar = X.(c_var{varInd});

selMat = mVar(mVar(:,SITE) == area & mVar(:,COND) == condition,:);

stimInd = ismember(selMat(:,POL),stimulation);
selMat = selMat(stimInd,:);

%Separate left and right initial hand conditions
selMat(:,REACH) = selMat(:,TARGET) - selMat(:,IHAND);

%For PMd remove KP
if (area == 2)
    selMat(selMat(:,SUB) == 7,:) = []; 
end

subject = unique(selMat(:,SUB)); 
%Draw a regression line between variable of interest and score of interest
if (gType == 1)
    
    for s = 1:length(unique(selMat(:,SUB)))
        
        %Select subject matrix
        subMat = selMat(selMat(:,SUB) == subject(s),:);
        
        %Compute regression on target 
        targMat = subMat(:,1:2); 
        
        TX = [ones(size(targMat,1),1) subMat(:,2)];
        b = regress(subMat(:,1),TX); 
        
        slope(s) = b(2); 
        
        
    end
    
    %Temporary due to automation 
    mu = slope; 
    
elseif (gType == 2)
    for s = 1:length(unique(selMat(:,SUB)))
        
        %Select subject matrix
        subMat = selMat(selMat(:,SUB) == subject(s),:);
        
        %Select left and right reach directions
        LReach = subMat(subMat(:,REACH) < 0,[1 8]);
        RReach = subMat(subMat(:,REACH) > 0,[1 8]);
        
        %Now construct a regression matrix
        LX = [ones(size(LReach,1),1) LReach(:,2)];
        RX = [ones(size(RReach,1),1) RReach(:,2)];
        
        %Perform regressions
        bL = regress(LReach(:,1),LX);
        bR = regress(RReach(:,1),RX);
        
        %Compute slope differences
        slope(s) = bL(2) - bR(2);
        
        %Compute means
        muL = mean(LReach(:,1));
        muR = mean(RReach(:,1));
        mu(s) = muL - muR;
    end
    slope = 0; 
elseif gType == 3
    %First select left IHP only
    selMat = selMat(selMat(:,IHAND) < 0,:);
    
    for s = 1:length(unique(selMat(:,SUB)))
        
        %Select subject matrix
        subMat = selMat(selMat(:,SUB) == subject(s),:);
        
        %Select left and right reach directions
        LHemi = subMat(subMat(:,TARGET) < 0,1);
        RHemi = subMat(subMat(:,TARGET) > 0,1);
        
        %Now construct a regression matrix
        Hemi = mean(LHemi) - mean(RHemi); 
        mu(s) = Hemi; 
    end
    slope = 0;
elseif (gType == 4) 
     for s = 1:length(unique(selMat(:,SUB)))
        
        %Select subject matrix
        subMat = selMat(selMat(:,SUB) == subject(s),:);
        
        %Select left and right reach directions
        LHand = subMat(subMat(:,IHAND) < 0,1);
        RHand = subMat(subMat(:,IHAND) > 0,1);
        
        %Now construct a regression matrix
        Hand = mean(LHand) - mean(RHand); 
        mu(s) = Hand; 
     end
    slope = 0;
elseif (gType == 5)
    for s = 1:length(unique(selMat(:,SUB)))
        %Select left IHP AND left Target
        subMat = selMat(selMat(:,SUB) == subject(s),:); 
        
        L = subMat(subMat(:,IHAND) < 0 & subMat(:,TARGET) < 0,1); 
        R = subMat(subMat(:,IHAND) > 0 & subMat(:,TARGET) > 0,1); 
        
        Hand = mean(L) - mean(R); 
        mu(s) = Hand; 
        slope = 0; 
        
    end
elseif (gType == 6)
    
    for s = 1:length(unique(selMat(:,SUB)))
        subMat = selMat(selMat(:,SUB) == subject(s),:); 
        
        %Grab left and right ihps
        L = subMat(subMat(:,IHAND) < 0,:); 
        R = subMat(subMat(:,IHAND) > 0,:); 
        
        LX = [ones(size(L,1),1) L(:,2)]; 
        RX = [ones(size(R,1),1) R(:,2)]; 
        
        bL = regress(L(:,1),LX); 
        bR = regress(R(:,1),RX); 
        
        slope(s) = bL(2) - bR(2); 
    end
    mu = slope; 
end
end





