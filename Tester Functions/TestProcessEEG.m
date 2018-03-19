function TestProcessEEG(E, T, alignType)
%PROCESSEEG Summary of this function goes here
%   Requires E Struct as an input. Processes information and outputs EEG
%   object
% Align types (0 or empty - target onset, 1 - movement onset)
%% Argument housekeeping

if nargin < 3
    alignType = 0;
end

%% Initialize Variables
a_str = {'Target Onset' 'Movement Onset'};
xl = [-0.5 -0.8];
X_MIN = xl(alignType + 1);
X_MAX = 1.0;
FREQ = 512;

%% Initialize EEG Parameters
subInitials = {E.subInitials};
subChar = cellfun(@(x) char(x), subInitials, 'un', 0);
disp(['Subject Initials: ' strjoin(subChar)]);

%Compute median reaction time in real time?
if alignType
    mu_rxn_samp = round(nanmedian(T{:,end})*FREQ/1000);
else
    mu_rxn_samp = 0;
end

%Get 90th percentile of movements
RT_col = T{:,end};
RT_col(isnan(RT_col)) = [];
q = [0.05 0.95];
nineP = round(quantile(RT_col,q)/1000*FREQ);

%% Perform preprocessing on data
filtT = [];
for s = 1:length(E)
    
    disp(s);
    S = E(s);
    TestPreProcessEEG(S, alignType, FREQ, mu_rxn_samp);
    
end

%Store all relevant information into EEG object
disp('Finished processing EEG.')

end

