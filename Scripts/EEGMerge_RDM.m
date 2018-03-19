%%Organization script for automated loading into EEG program.

%Use this program prior to EEG analysis on DMAT files with filename type:
%IIASCCBD.mat.mat
% Where:
% I is intials
% A is area of stimulation
% S is stimulation type
% C is condition type
% B is block number
%Written by - Jerry Jeyachandra 2016
%Output Structure - %S--> subjects --> A/C --> Pre/Post --> D Matrices/EEG Data

%TO-DO: Nothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get Directory
clear;
p = uigetdir('Select location of DMats...');
cd(p);

F = dir;

%Remove filesystem rows
F([F.isdir] == 1) = [];

%% Create all uppercase field to eliminate issue of case-sensitivity in many of MATLAB's functions
upperCell = upper({F.name});
upperStruct = cell2struct(upperCell,'uName',1);
[F(:).uName] = deal(upperStruct.uName);

%%  Initialize variables
%Indices based on standardized file name schema, change here if you have a
%different schema - rest of script should work fine.
ind_stim = 3;
ind_initial = 1:2;
ind_cond = 5:6;
ind_block = 7;
f_stim = {'Anodal' 'Cathodal'};
f_cond = {'Pre' 'Post' 'Stim'};
E = struct;
merge_T = [];

sp = uigetdir('Select location of saved MAT file...');
time = clock;
strFile = ['[' num2str(time(1)) '-' num2str(time(2)) '-' num2str(time(3)) '-' num2str(time(4)) '-' num2str(time(5)) ']' ];
save([sp '\EEGStructRDM' strFile], 'E', '-v7.3');



%% Sort Struct (just in case filesystem sorts weirdly...), MAC does weird stuff
FCell = struct2cell(F);
FFields = fields(F);
szF = size(FCell);
FCell = reshape(FCell, szF(1), []);
FCell = FCell';
FCell = sortrows(FCell,1);
FCell = reshape(FCell',szF);
FSorted = cell2struct(FCell,FFields,1);
cell_F = struct2cell(FSorted);
fNameC = {FSorted.uName};

%% Get all subject initials, stimulation types and conditions
%Retrieve subject initials
sub_initials = unique(cellfun(@(x) x(ind_initial), cell_F(end,:), 'un', 0))';

%Check that cond_type actually does pre then post in the indices, flip if
%not! Dependent on file naming convention!
cond_type = sortrows(unique(cellfun(@(x) x(ind_cond), cell_F(end,:), 'un', 0)))';
%(1*)
cond_type(end) = [];

%Check that anodal precedes cathodal!
stim_type = sortrows(unique(cellfun(@(x) x(ind_stim), cell_F(end,:), 'un', 0)))';

%Data structure for storing conditional indexing information
T = table;
T_test = [];
delCount = 0; %To account for removals due to early start

%Check for weird non-matching zeros to RT

zeroTest = [];
%% Organize struct
for s = 1:length(sub_initials)
    
    %Grab all files with subject initials
    subInd = ~cellfun('isempty',strfind(cellfun(@(x) x(ind_initial), fNameC, 'un', 0),sub_initials{s}));
    
    tMovementList = [];
    for stim = 1:length(stim_type)
        
        %Divide into anodal/cathodal conditions
        stimInd = ~cellfun('isempty',strfind(cellfun(@(x) x(ind_stim), fNameC, 'un', 0),stim_type{stim}));
        
        for cond = 1:length(cond_type)
            
            %Divide into pre/post conditions (since no STIM EEG)
            condInd = ~cellfun('isempty', strfind(cellfun(@(x) x(ind_cond), fNameC, 'un', 0), cond_type{cond}));
            
            %Logical containing struct indices pertaining to file with
            %specific sub,stim,cond.
            all_Ind = subInd & stimInd & condInd;
            cur_Struct = FSorted(all_Ind);
            
            %Initialize transient variables
            merge_S = [];
            merge_EEG = [];
            merge_END = [];
            merge_T = [];
            counter = 1;
            
            %Load files, merge DMATs, EEG and organize into struct
            
            for m = 1:sum(all_Ind)
                curFile = load(cur_Struct(m).name);
                
                
                %If missing EEG data in DMAT display filename for debugging
                if (~isfield(curFile.D{1},'eegData'))
                    disp(cur_Struct(m).name);
                end
                
                %Merge if EEG is complete or if trigger is valid 
                if (isfield(curFile.D{1},'eegData') && ~isnan(mean(curFile.D{1}.eegData(:,1))) && ~sum(find(curFile.D{1}.eegData(:,20) == 255)))
                    
                    %Merge EEG and assign block number
                    merge_EEG = [merge_EEG; curFile.D{1}.eegData(:,[1:6 20]) m*ones(size( curFile.D{1}.eegData,1),1)];
                    X = [curFile.D{2:end}];
                    
                    
                    
                    %Grab relevant events for EEG alignment
                    tFix = [curFile.D{1}.tFixation [X(:).tFixation]];
                    tMotOn = [curFile.D{1}.tMotOn [X(:).tMotOn]];
                    tMotEnd = [curFile.D{1}.tMotEnd [X(:).tMotEnd]];
                    tSPon = [curFile.D{1}.SPon [X(:).SPon]];
                    
                    %Short fix for misaligned EEG issues... look into files
                    %later...
                    Q = tMotEnd - tFix;
                    tMotOn(Q > 4000) = NaN;
                    tMotEnd(Q > 4000) = NaN;
                    
                    %Cumulative merged values across blocks
                    merge_T(:,:,counter) = [tFix' tMotOn' tSPon' + tFix'];
                    
                    if (size(merge_T,1) < 100) 
                        merge_T(size(merge_T,1):100,:,counter) = nan; 
                    end
                    
                    %Get saccade direction
                    eyeDec = [curFile.D{1}.eyedec [X(:).eyedec]];
                    tarDir = [curFile.D{1}.tdir [X(:).tdir]];
                    corrDec = [curFile.D{1}.correct [X(:).correct]];
                    
                    %Convert into saccade indexing
                    lrEye = (eyeDec == 180);
                    lrTarg = (tarDir == 180);
                    
                    %Indices for LR targ/saccade same direction, this
                    %implies correct though
                    L = (lrEye & lrTarg);
                    R = (~lrEye & ~lrTarg);
                    
                    %Convert into grouping indexers (L right vs wrong, R
                    %right vs wrong)
                    %Correct = 1, Wrong = 2
                    lCorr = 2 - (2.*R) - L;
                    rCorr = 2 - (2.*L) - R;
                    
                    %L vs R comparison (correct only, left = 1, right = 2)
                    LvR = L + 2.*R;
                    
                    %Correct vs Wrong trials
                    CvW = 2 - corrDec;
                    
                    %Convert matrix to table form (unidentified headers)
                    %Default set up:
                    %Subject, Stim, Cond, Block, Trial, Inds
                    numTrials = ones(size(tSPon',1),1); %Any variable size can be used...
                    T_add = table(s.*numTrials,stim.*numTrials,cond.*numTrials,...
                        m.*numTrials, [1:length(numTrials)]', lCorr', rCorr', LvR', CvW', tSPon');
                    
                   
                    
                    %Removal of bad trials set to NaN, EEG PreProcessing
                    %will take care of NaNs in the event timing matrix
                    %(merge_T which finally forms T matrix)
                    trialMark = [curFile.D{1}.good [X(:).good]];
                    badInd = trialMark == 1;
                    badSaccade = tSPon == 0;
                    merge_T(badInd | badSaccade,:,counter) = NaN;
                    
                    
                    %Search for NaNs applied and from actual recording
                    %errors
                    %Table needs NaNs removed here, since it is not thrown
                    %into the EEG preprocessing script
                    T_add(isnan(merge_T(1:height(T_add),1,counter)),:) = [];
                    
                    %Short fix for EEG start after trial begins
                    %Really should just set it to NaN, but this works
                    if (curFile.D{1}.eegData(1,20) == 3)
                        merge_T(1,:,counter) = 0;
                        T_add(1,:) = [];
                        delCount = delCount + 1;
                    end
                    
                    counter = counter + 1;
                    
                    T = vertcat(T,T_add);
                    T_test = [T_test; merge_T(:,:,counter-1)];
                    
                    %Routine 1: Check SPon, for RDM only
                    x = merge_T(:,end,counter-1) - merge_T(:,1,counter-1);
                    zeroTest = [zeroTest; x];
                    
                end
            end
            
            if (any(merge_T))
                %Construct struct containing organized info.
                E(s).(f_stim{stim}).(f_cond{cond}).EEG = merge_EEG;
                E(s).(f_stim{stim}).(f_cond{cond}).T = merge_T;
                
                %Test, build on current set. 
                %If you get an error here that means that the T table for
                %indexing groupings and T matrix for storing event timings
                %doesn't match. 
                %Fix this before running analysis!
                T_test(isnan(T_test(:,1)),:) = [];
                if (height(T)+delCount ~= size(T_test,1))
                    error('Critical Error: Non-matching sets!');
                end
            end
        end
        
    end
    %% Assign subject initial to each row of struct array
    E(s).subInitials = sub_initials(s);
    save([sp '\EEGStructRDM' strFile],'E', '-append');
end
%Stick on indexing table
%Filter out the last trial - this is due to EEGAlignEvent_v2 where last
%trials are skipped due to high variability in stopping time, if you stop
%too early after trial is over then you get an error. 
t50rows = T.Var5 == 100;
T(t50rows,:) = [];
save([sp '\EEGStructRDM' strFile], 'T','-append');
