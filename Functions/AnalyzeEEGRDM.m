function [ TFRhandler, dSTFR, Apre, Cpre, Apost, Cpost] = AnalyzeEEGRDM
% EEG Analysis written by Gunnar Blohm modified for Random Dot Motion Task
% USE WITH EEGMergeRDM OUTPUT
% Written by Jerry Jeyachandra - 2016.
%TO-DO: Add left/right comparisons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data if not already loaded (save time if already run previously)
if (~exist('E', 'var'))
    [EEGFile, pe] = uigetfile(cd, 'Select EEGStruct File...');
    load([pe EEGFile]);
end

clearvars -except E;

subInitials = {E.subInitials};
subChar = cellfun(@(x) char(x), subInitials, 'un', 0);
disp(['Subject Initials: ' strjoin(subChar)]);

%% Get Alignment type (to Target Onset or to Arm Movement Start)

validAInput = false;
while (~validAInput)
    disp('Target Onset - 1 \n Saccade Onset - 2');
    alignInput = input('Enter event to be aligned, enter ! to quit: ', 's');
    
    if (strcmpi(alignInput,'!'))
        error('Program Aborted');
    else
        alignBool = str2double(alignInput) - 1;
        validAInput = true;
    end
end


a_Str = {'Target Onset' 'Saccade Onset'};
alignStr = a_Str{str2double(alignInput)};
xl = [-0.5 -1.65];
X_MIN = xl(alignBool + 1);
X_MAX = 1.0;

%% Initialize variables
FIX_TIME = 516; %RDM experiment
Fs = 512;
ITI = 551;
Ts = 1/Fs;
lab = {'Target Onset', 'Saccade Onset'};
ApreTFR = [];
CpreTFR = [];
ApostTFR = [];
CpostTFR = [];
c_range = [1 50];

%% Get Median Latency across all conditions, stimulation types etc... 
Y = []; 
X = struct2cell(E);
X = X(1:2,:,:);
X = squeeze(X);
X = X(:);
X = cell2mat(X);
X = struct2cell(X);
X = cell2mat(X);
X = struct2cell(X);
X = X(2,:,:);
X = squeeze(X);
X = X(:);
for i = 1:length(X)
    Y = cat(3,Y,X{i});  
end
Z = Y(:,3,:); 
Z = Z(:); 
mu_rxn = median(Z); 
mu_rxn_samp = mu_rxn*Fs/1000;

%% Perform TFR analysis for all STIM and COND
for sub = 1:length(E);
    
    %Get current set of fields and struct
    S = E(sub);
    s_fields = fields(S);
    c_fields = fields(S.Anodal);
    
    for stim = 1:length(s_fields) - 1
        for cond = 1:length(c_fields)
            
            %Get number of blocks in EEG file
            validBlock = unique(S.(s_fields{stim}).(c_fields{cond}).EEG(:,end));
            for k = 1:5
                eeg{k} = [];
            end
            
            lrInd = [];
            if (isempty(lrInd))
                lrVec = 1; 
            else
                lrVec = 1:2; 
            end
            for b = 1:length(validBlock)
                
                %%Grab current block of EEG data
                cS = S.(s_fields{stim}).(c_fields{cond});
                dat = cS.EEG;
                data = dat(dat(:,end) == validBlock(b),:);
                data(:,2:6) = data(:,2:6)./sqrt(1000);
                
                %% Find Trials
                ind = find(data(:,7) == 3);
                trig = group(ind, 1, 10); %Default crit is 2
                
                %Grouping check
                if (length(trig) ~= 200)
                    disp(['Problematic grouping!!! \n b = ' num2str(b) ' \n c = ' num2str(cond) ' \n s = ' num2str(stim)]);
                    pause;
                end
                
                Teeg = 1000*(data(end,1) - data(1,1))/length(data(:,1));
                
                events = [];
                
                
                events(:,1) = cS.T(:,1,b)./1000; %Motion Start
                events(:,2) = cS.T(:,2,b)./1000; %Motion End
                events(:,3) = cS.T(:,3,b)./1000; %Saccade onset
                
                                
                %Remove first trial if trial began prior to EEG
                if (cS.T(1,:,validBlock(b)) == 0)
                    trig(1:2) = [];
                    events(1,:) = [];
                end
                
                switch(alignBool)
                    case 0
                        %Target onset - 516ms after EEG recordings
                        alignEvent = ones(size(events,1),1) .* FIX_TIME;
                        alignEvent_samp = round(alignEvent*Fs/1000);
                        alignBaseOffset = 0;
                    case 1
                        %
                        alignEvent = events(:,3);
                        alignEvent_samp = round(alignEvent*Fs);
                        alignBaseOffset = round(mu_rxn_samp);
                end
                
                %% Filtering and Pre-processing
                %Apply pre-processing separately to each block
                dataf = [];
                zdataf = [];
                probch = [];
                
                %Butterworth Filter Transformation Function (6th order)
                [nb, db] = butter(3, 2*Teeg./1000*[1 50]);
                [nn, dn] = butter(3, 2*Teeg./1000*[55 65], 'stop');
                
                %Filter
                for i = 1:5
                    
                    dataf(:,i) = data(:,i+1);
                    %Apply notch filter 55-65Hz exclusion
                    dataf(:,i) = filtfilt(nn,dn,data(:,i+1));
                    %Apply bandpass filter 1-50Hz inclusion
                    dataf(:,i) = filtfilt(nb,db,dataf(:,i));
                    %Normalize data (assuming normal dist?)
                    dataf(:,i) = (dataf(:,i) - mean(dataf(:,i)))./std(dataf(:,i));
                    %Compute zscore of consecutive differences
                    zdataf = abs(zscore(diff(data(:,i+1))));
                    %Find problem channels if exists, if sum of zscores exceeds 7.5
                    probch = find(sum(zdataf > 7.5) > 0);
                    
                end
                for i = 1:length(probch)
                    disp(['Channel ' num2str(probch(i)-1) ' is problematic']);
                end
                                
                %% Align with Time of Interest (task dependent)
                for k = 1:5
                    for i = 1:length(trig)/2
                        
                        %Grab trigger blocks
                        curTrialInd = ind(trig(2*i-1):trig(2*i));
                        
                        %Updated Code for Movement Onset capability
                        arr = curTrialInd(alignEvent_samp(i))-(250 + alignBaseOffset):curTrialInd(alignEvent_samp(i)) + 500;
                        
                        Beeg{k}(:,i) = dataf(arr,k);
                        
                    end
                end
                
                for i = 1:length(Beeg)
                    eeg{i} = [eeg{i} Beeg{i}];
                end
                
            end
            
            
            %% TFR Analysis
            %Laplacian Time Series
            EEG = eeg{1} + 0.25*(eeg{2}+eeg{3}+eeg{4}+eeg{5});
            
            
            
            %Separate into left and right reach trials
            for lr = lrVec
                
                %Grab left/right reach vector EEG set.
                if (isempty(lrInd))
                    dirEEG = EEG; 
                else
                    dirEEG = EEG(:,lrInd == lr);
                end
                
                
                Freq = 1:45; % frequency range
                Time = data(arr,1)-data(arr(251 + alignBaseOffset),1); % time range
                
                % use makeTFR2 for baseline variance-normalized TFR = z-scores
                tfr{1} = makeTFR2(dirEEG, Freq, Time, 1000/Teeg, 7, data([arr(1) arr(251)],1)-data(arr(251),1));
                TFR(:,:, lr, cond) = tfr{1};
                
                [x, y] = meshgrid(Time, Freq);

            end
        end
        
        %% Compute TFR Change
        dTFR(:,:,:,stim) = TFR(:,:,:,2) - TFR(:,:,:,1);        
        
        %Get Pre/Post TFR
        subpreTFR(:,:,:,stim) = TFR(:,:,:,1);
        subpostTFR(:,:,:,stim) = TFR(:,:,:,2);
        
    end
    
    %Concatenate subject TFRs
    ApreTFR = cat(4, ApreTFR, subpreTFR(:,:,:,1));
    CpreTFR = cat(4, CpreTFR, subpreTFR(:,:,:,2));
    ApostTFR = cat(4,ApostTFR,subpostTFR(:,:,:,1));
    CpostTFR = cat(4,CpostTFR,subpostTFR(:,:,:,2));
    
    %Plot Difference TFR for each subject
    dSTFR(:,:,:,sub) = dTFR(:,:,:,1) - dTFR(:,:,:,2);
    
end

Apre = ApreTFR; 
Cpre = CpreTFR; 
Apost = ApostTFR; 
Cpost = CpostTFR; 


%Initialize TFR plotting handler object 
TFRhandler = TFRPlotHandler(x,y,X_MIN,X_MAX,Freq,Time(1),Time(251)); 

end

