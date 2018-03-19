classdef EEGMatrixHandler
    %EEGMATRIXHANDLER handles all EEG data analysis in one file in order to
    %ease extension, accessibilty and transferrability of software.
    
    %When using for the first time set up SetMode method to match your
    %grouping categories first! Give your experiment it's own mode number
    %1 -- RDM
    %2 -- Arm Reach
    %3 -- Free Choice    
    
    properties
        data; %rows - stimulation %columns - condition
        T; 
        channels;
        a_type; %alignment type
        Fs; %sampling frequency
        Time; %time vector (for plotting)
        Freq; %Frequencies for TFR
        surfX; %Plotting surface X
        surfY; %Plotting surface Y
        NUM_COND = 11; %columns of LRVEC indexer - 10 without RT
        LEFT = 1;
        RIGHT = 2;
        mode; %1 - Saccade, 2 - Arm Reach
        savePath; %Specified savePath
        baseOffset; %mean reaction time
        medianTime; 
        nP; % 95% movement range
    end
    
    methods
        
        %% Constructor and initialization
        function obj = EEGMatrixHandler(data, channels, alignType, Fs, baseOffset,nineP,T)
            
            obj.data = data;
            obj.a_type = alignType;
            obj.Fs = Fs;
            obj.Time = ((1:size(data{1,1,1},2))*(1/Fs)) - (251+baseOffset)*(1/Fs);
            obj.Freq = 1:45;
            [obj.surfX, obj.surfY] = meshgrid(obj.Time,obj.Freq);
            obj.channels = channels;
            obj.savePath = cd;
            obj.baseOffset = baseOffset;
            
            
            %Fix for nP exceeding Time vector indices
            %If less than time use nineP, if longer than time, use obj.Time
            %(subracted by 251, to remove baseline inclusion)
            nineP(2) = nineP(2) * (nineP(2) < length(obj.Time)) ...
                + (length(obj.Time)-251) * (nineP(2) > length(obj.Time));
            
            obj.nP = nineP;
            obj.mode = 2; %CHANGE BACK  
            
            obj.T = T;
            obj.T.Properties.VariableNames = obj.SetMode(obj.mode);
            obj.medianTime = median(obj.T.RT); 
            
        end
        
        %Set up the table variables based on the experiment
        %RT is treated as a special case in the analysis, explicitly use RT
        %as a name if you'd like to perform RT quantile separation and
        %comparisons 
        function varNames = SetMode(~,EMode)
            switch (EMode)
                case 1 %RDM
                    varNames = {'Subject', 'Stim', 'Cond', 'Block' ...
                        ,'Trial','L','R','LR','CW','RT'};
                case 2 %ARMREACH
                    varNames = {'Subject', 'Stim', 'Cond', 'Block' ...
                        ,'Trial', 'TIHP', 'targIHP', 'reachVec', 'EHemifield', 'LHemifield', 'RHemifield', ...
                        'Reach', 'Target', 'None', 'IHP', 'RT'};
                case 3 %FC
                    varNames = {'Subject', 'Stim', 'Cond', 'Block' ...
                        ,'Trial'};
            end
            
        end
        
        function obj = RemoveLearning(obj)
            
            obj.T = RemoveT_LR(obj.T); 
            
        end
        
        %% General Data Accessing
        
        %Returns specified matrix irrespective of grouping
        function selMat = SelectMatrix(obj,subject,condition,stimulation, bLaplace)
            
            selMat = vertcat(obj.data{stimulation,condition,subject});
            if (bLaplace)
                selMat = selMat - vertcat(obj.channels{stimulation,condition,subject});
            end
            
            selMat(:,end-obj.NUM_COND:end) = [];
            
        end
        
        %Returns grouped matrices
        function [lSelMat, rSelMat] = SelectGroupMatrix(obj,subject,condition,stimulation, bLaplace, group)
            
            selMat = vertcat(obj.data{stimulation,condition,subject});
            
            %Spatial Laplace filter
            if (bLaplace)
                selMat = selMat - vertcat(obj.channels{stimulation,condition,subject});
            end
            
            %Select subject, condition and stimulation of table
            subInd = ismember(obj.T.Subject,subject);
            stimInd = ismember(obj.T.Stim,stimulation);
            condInd = ismember(obj.T.Cond,condition);
            selInd = subInd & stimInd & condInd;
            
            %Define new table containing values for specific Sub/Cond/Stim
            selT = obj.T(selInd,:);
            
            % BUGTESTING 
%             if (height(selT) ~= size(selMat,1))
%                 disp([subject, stimulation, condition]); 
%             end
                        
            if (strcmpi(group,'None'))
                gLInd = logical(ones(height(selT),1));
                lSelMat = selMat(gLInd,:); 
                rSelMat = [];
            else
                
                %Left/Right indexers
                gLInd = selT.(group) == obj.LEFT;
                gRInd = selT.(group) == obj.RIGHT;
                
                %Grab left and right matrices and return
                lSelMat = selMat(gLInd,:);
                rSelMat = selMat(gRInd,:);
            end
            
            
            
            
        end
        
        %Returns binned structure for a given grouping parameter (RT for
        %now, expand if needed)
        function binMat = SelBinMat(obj,subject,condition,stimulation,binSiz,bLaplace)
            
            selMat = vertcat(obj.data{stimulation,condition,subject});
            selMat = selMat - vertcat(obj.channels{stimulation,condition,subject});
            %Select subject, condition and stimulation of table
            subInd = ismember(obj.T.Subject,subject);
            stimInd = ismember(obj.T.Stim,stimulation);
            condInd = ismember(obj.T.Cond,condition);
            selInd = subInd & stimInd & condInd;
            
            %Define a new table with specific Sub/Cond/Stim
            selT = obj.T(selInd,:);
            
            %Break down into trial % bins, sizes should be roughly equal
            pBin = linspace(0,1,binSiz+1);
            
            %Sort table based on reaction times, sort indices for selmat in
            %the same way
            [sortT,indT] = sortrows(selT,'RT');
            sortMat = selMat(indT,:);
            
            %Bin trials according to parameter set
            numTrials = size(sortMat,1);
            for i = 1:length(pBin)-1
                binMat{i} = sortMat(ceil(pBin(i)*numTrials)+1 : ceil(pBin(i+1)*numTrials),:);
            end
            
        end
        
        function varNames = GetGroupNames(obj)
            varNames = obj.T.Properties.VariableNames(6:end); 
        end
        
        %% TFR
        
        %Compute TFR
        function TFR = ComputeTFR(obj,subject,condition,stimulation, varargin)
            
            validTypes = {'Pseudo-Z', 'Power'};
            validLaplace = {'None', 'Laplacian'};
            if (obj.mode == 1)
                validGroup = {'None' 'Saccade' 'Target'};
            else
                validGroup = {'IHP' 'None' 'Target' 'Reach' 'RHemifield' 'LHemifield' 'EHemifield' 'MedianSplit'};
            end
            
            Results = parseInput(obj.mode,varargin{:});
            
            TFRtype = find(strcmpi(validTypes,Results.TFRtype));
            bLaplace = find(strcmpi(validLaplace,Results.SpatialFilter)) - 1;
            
            if (~strcmpi(Results.Group,'None'))
                [lMat, rMat] = obj.SelectGroupMatrix(subject,condition,stimulation,bLaplace,Results.Group);
                lTFRMat = MakeTFRP(lMat(:,1:end)',obj.Freq,obj.Time(1:end-obj.NUM_COND),obj.Fs,7,[obj.Time(1) obj.Time(251)], TFRtype);
                lTFR = TFRObject(lTFRMat, TFRtype, subject, condition, stimulation, bLaplace,obj.LEFT);
                
                rTFRMat = MakeTFRP(rMat(:,1:end)',obj.Freq,obj.Time(1:end-obj.NUM_COND),obj.Fs,7,[obj.Time(1) obj.Time(251)], TFRtype);
                rTFR = TFRObject(rTFRMat, TFRtype, subject, condition, stimulation, bLaplace,obj.RIGHT);
                TFR = [lTFR, rTFR];
                
            else
                curMat = obj.SelectMatrix(subject,condition,stimulation, bLaplace);
                TFRMat = MakeTFRP(curMat(:,1:end)',obj.Freq,obj.Time(1:end-obj.NUM_COND),obj.Fs,7,[obj.Time(1) obj.Time(251)], TFRtype);
                TFR = TFRObject(TFRMat, TFRtype, subject, condition, stimulation, bLaplace);
            end
            
        end
        
        %Display TFR
        function hFig = DispTFR(obj, TFR, varargin)
            
            
            p = inputParser;
            
            defaultRange = [TFR.GetRange*-1 TFR.GetRange];
            checkRange = @(x) isnumeric(x);
            addParameter(p,'CAxis',defaultRange,checkRange);
            
            p.KeepUnmatched = true;
            parse(p,varargin{:});
            
            c_range = p.Results.CAxis;
            
            a_Str = {'Target Onset' 'Movement Onset'};
            t_str = {'Pseudo-Z' 'Power'};
            d_str = {'Left' 'Right'};
            
            alignStr = a_Str{obj.a_type + 1};
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            figure;
            surf(obj.surfX,obj.surfY, TFR.GetTFR);
            shading interp
            view([0 0 1])
            
            
            caxis([c_range(1) c_range(2)])
            
            
            title([alignStr ' Direction ' d_str{TFR.GetDirection}]);
            cc = colorbar;
            ylabel(cc, t_str{TFR.GetType})
            colormap('Jet');
            hold on;
            
            xlim([X_MIN X_MAX])
            
            plot3([obj.Time(1) obj.Time(1)], [0 max(obj.Freq)], [2 2], 'k' );
            plot3([obj.Time(251 + obj.baseOffset) obj.Time(251 + obj.baseOffset)], [0 max(obj.Freq)], [2 2], 'k' );
            ylabel('Frequency (Hz)')
            xlabel(['Time from ' a_Str ' (ms)']);
            
            if (obj.a_type == 0)
                %Include patch for 90% movement
                timeNP = obj.Time(obj.nP(1)+251:obj.nP(2)+251);
                plot3([timeNP(1) timeNP(1)],[0,45],[2 2],'k:'); 
            end
            
            
            hFig = gcf;
            
            
        end
        
        function cFig = OverlayRS(obj,rsMat, thres) 
            
            if nargin < 3
                thres = 0.05 
            end
            
            rsBinary = rsMat < thres; 
            
            %Get significant clusters 
            [B, L] = bwboundaries(rsBinary); 
            for i = 1 : max(length(B))
                
                boundary = B{i}; 
                Z = 2.*ones(size(boundary,1),1); 
                
                plot3(obj.Time(boundary(:,2)),boundary(:,1),Z,'k-','LineWidth',1); hold on;
                
            end            
            
        end
        
        %Save TFR from Object (for post TFR computations)
        function SaveCTFR(obj,TFR,saveStr, varargin)
            if nargin < 3
                saveStr = '';
            end
            
            c_cond = {'Pre' 'Post' 'DC'};
            s_cond = {'A' 'C' 'DS'};
            
            hFig = obj.DispTFR(TFR, varargin{:});
            
            saveas(hFig, fullfile(obj.savePath,['TFR_s' num2str(TFR.GetSubject) ...
                '_t' s_cond{TFR.GetStimulation} '_c' c_cond{TFR.GetCondition} saveStr]),'fig');
            
            
        end
        
        %Subtracts TFR2 from TFR1
        function TFR = GetDifference(~, TFR1, TFR2)
            
            %Check for conditional subtraction or stimulation subtraction
            if (TFR1.GetStimulation ~= TFR2.GetStimulation)
                stimulation = 3;
            else
                stimulation = TFR1.GetStimulation;
            end
            
            if (TFR1.GetCondition ~= TFR2.GetCondition)
                condition = 3;
            else
                condition = TFR1.GetCondition;
            end
            
            TFR = TFRObject(TFR1.GetTFR-TFR2.GetTFR,TFR1.GetType,TFR1.GetSubject,condition,stimulation, TFR1.GetFilter);
        end
        
        %Vertically collapse array
        function TFRArray = GetArrayDifference(obj,TFRArrayIn)
            for i = 1:size(TFRArrayIn,2)
                TFRArray(i) = obj.GetDifference(TFRArrayIn(1,i),TFRArrayIn(2,i));
            end
        end
        
        %Save TFR
        function SaveTFR(obj,subject,condition,stimulation,varargin)
            
            p = inputParser;
            
            defaultType = 'Pseudo-Z';
            validTypes = {'Pseudo-Z', 'Power'};
            checkType = @(x) any(validatestring(x,validTypes));
            
            defaultLaplace = 'None';
            validLaplace = {'None', 'Laplacian'};
            checkLaplace = @(x) any(validatestring(x,validLaplace));
            
            addParameter(p, 'TFRtype', defaultType, checkType);
            addParameter(p, 'SpatialFilter', defaultLaplace, checkLaplace);
            
            p.KeepUnmatched = true;
            
            parse(p,varargin{:});
            
            TFR = obj.ComputeTFR(subject,condition,stimulation, p.Results.TFRtype, p.Results.SpatialFilter);
            
            hFig = obj.DispTFR(TFR);
            
            c_cond = {'Pre' 'Post'};
            s_cond = {'A' 'C'};
            saveas(hFig, fullfile(obj.savePath,['TFR_s' num2str(subject) ...
                '_t' s_cond{stimulation} '_c' c_cond{condition}]),'fig');
            
        end
        
        %Save all subjects TFR
        function SaveAllTFR(obj, varargin)
            
            p = inputParser;
            
            defaultType = 'Pseudo-Z';
            validTypes = {'Pseudo-Z', 'Power'};
            checkType = @(x) any(validatestring(x,validTypes));
            
            defaultLaplace = 'None';
            validLaplace = {'None', 'Laplacian'};
            checkLaplace = @(x) any(validatestring(x,validLaplace));
            
            addParameter(p, 'TFRtype', defaultType, checkType);
            addParameter(p, 'SpatialFilter', defaultLaplace, checkLaplace);
            
            p.KeepUnmatched = true;
            
            parse(p,varargin{:});
            
            c_cond = {'Pre' 'Post'};
            s_cond = {'A' 'C'};
            
            for s = 1:size(obj.data,3)
                for t = 1:size(obj.data,1)
                    for c = 1:size(obj.data,2)
                        TFR = obj.ComputeTFR(s,c,t,p.Results.TFRtype, p.Results.SpatialFilter);
                        hFig = obj.DispTFR(TFR);
                        
                        saveas(hFig, fullfile(obj.savePath,['TFR_s' num2str(s) ...
                            '_t' s_cond{t} '_c' c_cond{c}]),'fig');
                    end
                end
            end
            
        end
        
        %Create graph of all subjects EEG
        function saveAllEEG(obj,varargin)
            p = inputParser;
            
            defaultLaplace = 'None';
            validLaplace = {'None', 'Laplacian'};
            checkLaplace = @(x) any(validatestring(x,validLaplace));
            
            addParameter(p, 'SpatialFilter', defaultLaplace, checkLaplace);
            
            p.KeepUnmatched = true;
            
            parse(p,varargin{:});
            
            
            c_cond = {'Pre' 'Post'};
            s_cond = {'A' 'C'};
            
            for s = 1:size(obj.data,3)
                for t = 1:size(obj.data,1)
                    for c = 1:size(obj.data,2)
                        hFig = obj.dispEEG(s,c,t,p.Results.SpatialFilter);
                        saveas(hFig, fullfile(obj.savePath,['EEG_s' num2str(s) ...
                            '_t' s_cond{t} '_c' c_cond{c}]),'fig');
                    end
                end
            end
            
        end
        
        %Compute correlations
        function [S] = ComputeCorr(~, Obj_TFR1, Obj_TFR2)
            
            TFR1_all = [];
            TFR2_all = [];
            
            
            for s = 1:size(Obj_TFR1,2)
                TFR1_all = cat(3,TFR1_all,Obj_TFR1(s).GetTFR);
                TFR2_all = cat(3,TFR2_all,Obj_TFR2(s).GetTFR);
                
            end
            
            S = CorrTFR(TFR1_all,TFR2_all);
            
        end
        
        %Display correlational analysis results
        function [hFig, hBFig] = DispCorr(obj, S, boolThres)
            
            if nargin < 3
                boolThres = 0;
            end
            
            if (boolThres)
                S.r(S.P > 0.05) = NaN;
                S.beta(S.P > 0.05) = NaN;
            end
            
            a_Str = {'Target Onset' 'Movement Onset'};
            alignStr = a_Str{obj.a_type + 1};
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            figure;
            surf(obj.surfX,obj.surfY, S.r);
            shading interp
            view([0 0 1])
            
            caxis([-1 1])
            title(alignStr);
            cc = colorbar;
            ylabel(cc, 'Correlation Coefficient')
            colormap('Jet');
            hold on;
            
            xlim([X_MIN X_MAX])
            
            plot3([obj.Time(1) obj.Time(1)], [0 max(obj.Freq)], [2 2], 'k' );
            plot3([obj.Time(251) obj.Time(251)], [0 max(obj.Freq)], [2 2], 'k' );
            ylabel('Frequency (Hz)')
            xlabel(['Time from ' a_Str ' (ms)']);
            
            hFig = gcf;
            
            figure;
            surf(obj.surfX,obj.surfY, S.beta);
            shading interp
            view([0 0 1])
            
            title(alignStr);
            cc = colorbar;
            ylabel(cc, 'Beta Value')
            colormap('Jet');
            hold on;
            caxis([-5 5]);
            xlim([X_MIN X_MAX])
            
            plot3([obj.Time(1) obj.Time(1)], [0 max(obj.Freq)], [2 2], 'k' );
            plot3([obj.Time(251) obj.Time(251)], [0 max(obj.Freq)], [2 2], 'k' );
            
            xlabel(['Time from ' a_Str ' (ms)']);
            ylabel('Frequency (Hz)')
            
            hBFig = gcf;
            
            
            
        end
        
        %Display correlational analysis results
        function [hFig, hBFig] = DispCoeff(obj, coeffMat)
            
            a_Str = {'Target Onset' 'Movement Onset'};
            alignStr = a_Str{obj.a_type + 1};
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            figure;
            surf(obj.surfX,obj.surfY, coeffMat);
            shading interp
            view([0 0 1])
            
            caxis([-1 1])
            title(alignStr);
            cc = colorbar;
            ylabel(cc, 'Correlation Coefficient')
            colormap('Jet');
            hold on;
            
            xlim([X_MIN X_MAX])
            
            plot3([obj.Time(1) obj.Time(1)], [0 max(obj.Freq)], [2 2], 'k' );
            plot3([obj.Time(251) obj.Time(251)], [0 max(obj.Freq)], [2 2], 'k' );
            ylabel('Frequency (Hz)')
            xlabel(['Time from ' a_Str ' (ms)']);
        end
        
        %Get Mean of TFR Objects
        function Obj_TFR = TFRMean(~,TFRArray)
            TFR = zeros(size(TFRArray(1).GetTFR));
            
            for i = 1:length(TFRArray)
                TFR = TFR + TFRArray(i).GetTFR;
            end
            
            mTFR = TFR/length(TFRArray);
            
            Obj_TFR = TFRObject(mTFR,TFRArray(1).GetType, 0, TFRArray(1).GetCondition, ...
                TFRArray(1).GetStimulation, TFRArray(1).GetFilter, TFRArray(1).GetDirection);
            
        end
        
        %Get Variance of TFR Objects
        function Obj_TFR = TFRVar(obj,TFRArray)
            
            mTFR = obj.TFRMean(TFRArray);
            sqTFR = zeros(size(TFRArray(1).GetTFR));
            for i = 1:length(TFRArray)
                sqTFR = sqTFR + (TFRArray(i).GetTFR - mTFR.GetTFR).^2;
            end
            
            varTFR = sqTFR / length(TFRArray);
            
            Obj_TFR = TFRObject(varTFR,TFRArray(1).GetType, 0, TFRArray(1).GetCondition, ...
                TFRArray(1).GetStimulation, TFRArray(1).GetFilter, TFRArray(1).GetDirection);
            
        end
        
        %Compute for all subjects
        function TFRArray = ComputeAllTFR(obj, condition,stimulation, varargin)
            
            for i = 1:size(obj.data,3)
                TFRArray(:,i) = obj.ComputeTFR(i,condition,stimulation,varargin{:});
            end
            
        end
        
        %Perform a RankSum test using baseline period
        function RSMatrix = RankSum(~,TFRArray)
            
            baseSet = 251;
            
            %Get subject TFR arrays
            for i = 1:length(TFRArray)
                TFR(:,:,i) = TFRArray(i).GetTFR;
            end
            
            disp('Computing Ranksum...')
            
            allBaseTFR = TFR(:,1:baseSet,:);
            baseTFR = squeeze(mean(allBaseTFR,2));
            
            %Compute rank sum for each frequency
            for i = 1:size(TFR,1)
                
                basePixel = baseTFR(i,:);
                basePixel = basePixel(:);
                
                for j = 1:(size(TFR,2))
                    %Grab freq/timebin to be tested across all subjects
                    testPixel = TFR(i,j,:);
                    testPixel = testPixel(:);
                    
                    %Get the ranksum test of pixel
                    rPixel(i,j) = ranksum(testPixel,basePixel);
                end
            end
            
            RSMatrix = rPixel;
            
            
        end
        
        %Apply indexer to TFRs
        function TFRFilt = ApplyIndexFilt(~,TFRArray, IndFilter)
            
            for i = 1:length(TFRArray)
                TFRArray(i).TFR(~IndFilter) = NaN;
            end
            
            TFRFilt = TFRArray;
        end
        
        %Restricted window test
        function pTFR = TestTFR(obj,TFRArray,window)
            
            %If window not specified, perform test on all pixel units
            if (nargin < 2)
                window = 1:length(obj.Time);
            end
            
            for i = 1:length(TFRArray)
                TFR(:,:,i) = TFRArray(i).GetTFR;
            end
            
            %Create a building matrix
            pTFR = [];
            
            for i = window(1):window(end)
                %Get current time slot across all frequencies
                timeSlice = squeeze(TFR(:,i,:))';
                [~,sliceP] = ttest(timeSlice);
                pTFR = [pTFR sliceP'];
            end
            
            
            %Fill in NaNs in out of window values such that only values are
            %defined within the window specified
            nanWindow = nan(45,length(obj.Time));
            nanWindow(:,window) = pTFR;
            pTFR = nanWindow;
            
        end
        %% Phase Analysis
        
        %Compute and output a raw phase TFR object (depreciated) 
        function PLOC = ComputePhaseOld(obj,subject,condition,stimulation,varargin)
            
            validTypes = {'Pseudo-Z', 'Power'};
            validLaplace = {'None', 'Laplacian'};
            if (obj.mode == 1)
                validGroup = {'None' 'Saccade' 'Target'};
            else
                validGroup = {'IHP' 'None' 'Target' 'Reach' 'RHemifield' 'LHemifield' 'EHemifield' 'MedianSplit'};
            end
            
            Results = parseInput(obj.mode,varargin{:});
            
            TFRtype = find(strcmpi(validTypes,Results.TFRtype));
            bLaplace = find(strcmpi(validLaplace,Results.SpatialFilter)) - 1;
            groupType = find(strcmpi(validGroup,Results.Group));
            
            if (~strcmpi(Results.Group,'None'))
                [lMat, rMat] = obj.SelectGroupMatrix(subject,condition,stimulation,bLaplace,groupType);
                [~,lPLOCMat] = MakeTFRP(lMat(:,1:end)',obj.Freq,obj.Time(1:end-obj.NUM_COND),obj.Fs,7,[obj.Time(1) obj.Time(251)], TFRtype);
                lTFR = TFRObject(lPLOCMat, TFRtype, subject, condition, stimulation, bLaplace,obj.LEFT);
                
                [~,rPLOCMat] = MakeTFRP(rMat(:,1:end)',obj.Freq,obj.Time(1:end-obj.NUM_COND),obj.Fs,7,[obj.Time(1) obj.Time(251)], TFRtype);
                rTFR = TFRObject(rPLOCMat, TFRtype, subject, condition, stimulation, bLaplace,obj.RIGHT);
                PLOC = [lTFR, rTFR];
                
            else
                curMat = obj.SelectMatrix(subject,condition,stimulation, bLaplace);
                [~,PLOCMat] = MakeTFRP(curMat(:,1:end)',obj.Freq,obj.Time(1:end-obj.NUM_COND),obj.Fs,7,[obj.Time(1) obj.Time(251)], TFRtype);
                PLOC = TFRObject(PLOCMat, TFRtype, subject, condition, stimulation, bLaplace);
            end
            
        end
        
        %Compute and output phase trajectory TFR object 
        function PhaseArray = ComputePhase(obj,subject,condition,stimulation,varargin)
            
            p = inputParser;
            
            defaultGroup = 'Target';
            validGroup = obj.T.Properties.VariableNames(6:end);
            checkGroup = @(x) any(validatestring(x,validGroup));
            addParameter(p, 'Group', defaultGroup, checkGroup);
            
            p.KeepUnmatched = true;
            parse(p,varargin{:});
            
            specMat = obj.SpecifyMatrix(subject,condition,stimulation,p.Results.Group);
            
            %Conditionals (Using Pseudo-Z) 
            for i = 1:length(specMat);
                PhaseMat = PhaseTFR(specMat{i}',obj.Freq,obj.Time,...
                    obj.Fs,7,[obj.Time(1) obj.Time(251)]);

                %Apply appending algorithm 
                lin_PhaseMat = LinearizePhase(PhaseMat); 
                                
                PhaseArray(1,i) = TFRObject(lin_PhaseMat, 1, subject, condition, stimulation, 1);
            end
            
        end
        
        %Compute and output phase trajectories across all TFR objects
        function PhaseArray = ComputeAllPhase(obj,condition,stimulation,varargin) 
            
            for i = 1:size(obj.data,3)
                PhaseArray(:,i) = obj.ComputePhase(i,condition,stimulation,varargin{:});
            end
            
        end
        
        function AvgPhase = AveragePhase(~,PhaseArray) 
            
            for i = 1 : size(PhaseArray,2)
                
                %Extract TFR array 
                pTFR(:,:,i) = PhaseArray(i).GetTFR; 
                             
            end
            
            %Compute average mean for each element 
            for freq = 1 : size(pTFR,1)
                for samp = 1 : size(pTFR,2)
                    avg_vec = pTFR(freq,samp,:); 
                    mu_phase(freq,samp) = cirmean(avg_vec); 
                end
            end
            
            %Run the appending algorithm 
            lin_mu_phase = LinearizePhase(mu_phase); 
            
            %Output the phase object 
            AvgPhase = TFRObject(lin_mu_phase,1,0,PhaseArray(1).GetCondition, ...
                PhaseArray(1).GetStimulation,1); 
            
        end
        
        %Display phase shifts
        function DispPhase(obj,PhaseTFR)
            
            a_Str = {'Target Onset' 'Movement Onset'};
            t_str = {'Pseudo-Z' 'Power'};
            d_str = {'Left' 'Right'};
            
            alignStr = a_Str{obj.a_type + 1};
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;

            figure;
            surf(obj.surfX,obj.surfY,PhaseTFR.GetTFR);
            shading interp
            view([0 0 1])
                                    
            title([alignStr ' Direction ' d_str{PhaseTFR.GetDirection}]);
            hold on;
            
            xlim([X_MIN X_MAX])
            
            plot3([obj.Time(1) obj.Time(1)], [0 max(obj.Freq)], [2 2], 'k' );
            plot3([obj.Time(251 + obj.baseOffset) obj.Time(251 + obj.baseOffset)], [0 max(obj.Freq)], [2 2], 'k' );
            ylabel('Frequency (Hz)')
            xlabel(['Time from ' a_Str ' (ms)']);
            hFig = gcf;
            colorbar;
            colormap jet;
            
        end
        
        %Display inter-trial coherence
        function DispITC(obj,TFR)
            
            a_Str = {'Target Onset' 'Movement Onset'};
            t_str = {'Pseudo-Z' 'Power'};
            d_str = {'Left' 'Right'};
            
            alignStr = a_Str{obj.a_type + 1};
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            ITC = TFR.GetTFR;
            
            figure;
            surf(obj.surfX,obj.surfY, ITC);
            shading interp
            view([0 0 1])
            
            c_range = [0 0.5];
            caxis([c_range(1) c_range(2)])
            
            title([alignStr ' Direction ' d_str{TFR.GetDirection}]);
            cc = colorbar;
            ylabel(cc, 'ITC')
            colormap('Jet');
            hold on;
            
            xlim([X_MIN X_MAX])
            
            plot3([obj.Time(1) obj.Time(1)], [0 max(obj.Freq)], [2 2], 'k' );
            plot3([obj.Time(251 + obj.baseOffset) obj.Time(251 + obj.baseOffset)], [0 max(obj.Freq)], [2 2], 'k' );
            ylabel('Frequency (Hz)')
            xlabel(['Time from ' alignStr ' (ms)']);
            hFig = gcf;
            
        end
        
        %Compute inter-trial coherence
        function ITC = ComputeITC(~,TFR)
            
            phAng = TFR.GetTFR;
            euForm = zeros(size(phAng,1),size(phAng,2));
            for j = 1:size(phAng,3)
                
                %Compute exponent for each phase TFR and accumulate
                euForm = euForm + exp(1i * phAng(:,:,j));
                
            end
            ITC = TFR;
            ITC.TFR = abs(euForm)./size(phAng,3);
            
        end
        
        %Compute across binned RTs
        function [mITC,pMat] = ComputeRTITC(obj,subject,condition,stimulation,binSiz,varargin)
            
            result = parseInput(obj.mode,varargin{:});
            
            %Break down data into bins
            binMat = obj.SelBinMat(subject,condition,stimulation,binSiz,result.SpatialFilter);
            
            %Compute ITC given RT bins
            for i = 1:length(binMat)
                
                %Compute the phase angle progression for each trial
                [~,PLOCMat] = MakeTFRP(binMat{i}',obj.Freq,obj.Time,obj.Fs,7,[obj.Time(1) obj.Time(251)], 1);
                PLOC = TFRObject(PLOCMat, 1, subject, condition, stimulation, 0);
                
                %Convert this into an ITC plot
                ITC(i) = obj.ComputeITC(PLOC);
                %Save it into a test matrix
                pMat(:,:,:,i) = PLOC.GetTFR;
                
            end
            
            mITC = ITC;
            
        end
        
        %Perform binned RT computation across all subjects
        function oITC = ComputeAllRTITC(obj,condition,stimulation,binSiz,varargin)
            %Compute ITC given RT bins across all subjects
            %Return as a sub x bin matrix for subsequent computations
            %(regression, nope its circular, wouldn't that suggest bigger
            %magnitude?)
        end
        
        %% EEG
        
        %Display EEG (preprocessed signal from a subject)
        function hFig = dispEEG(obj, subject, condition, stimulation, varargin)
            
            results = parseInput(obj.mode, varargin{:});
            
            p = inputParser;
            
            defaultLaplace = 'None';
            validLaplace = {'None', 'Laplacian'};
            checkLaplace = @(x) any(validatestring(x,validLaplace));
            
            addParameter(p, 'SpatialFilter', defaultLaplace, checkLaplace);
            
            p.KeepUnmatched = true;
            
            parse(p,varargin{:});
            
            bLaplace = find(strcmpi(validLaplace,results.SpatialFilter)) - 1;
            
            meanMat = obj.getEEGMean(subject,condition,stimulation, bLaplace);
            varMat = obj.getEEGVar(subject,condition,stimulation, bLaplace);
            posVarMat = meanMat + varMat;
            negVarMat = meanMat - varMat;
            
            figure;
            plot(obj.Time,meanMat,'b', 'LineWidth', 1.5);
            hold on;
            plot([obj.Time(1) obj.Time(end)],[0 0], 'k--', 'LineWidth', 1); hold on;
            patch([obj.Time fliplr(obj.Time)], [posVarMat fliplr(negVarMat)], 'c', 'FaceAlpha', 0.5, 'EdgeColor', 'w');
            if (~bLaplace)
                y = ylim;
            else
                y = [-0.1 0.1];
            end
            plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5);
            hFig = gcf;
        end
        
        %Save a specific EEG
        function saveEEG(obj, subject, condition, stimulation, varargin)
            
            hFig = obj.dispEEG(subject, condition, stimulation, varargin{:});
            
            c_cond = {'Pre' 'Post'};
            s_cond = {'A' 'C'};
            saveas(hFig, fullfile(obj.savePath,['EEG_s' num2str(subject) ...
                '_t' s_cond{stimulation} '_c' c_cond{condition}]),'fig');
            
        end
        
        %Compute EEG mean
        function meanMat = getEEGMean(obj, subject, condition, stimulation, bLaplace)
            
            mat = obj.SelectMatrix(subject,condition,stimulation, bLaplace);
            meanMat = mean(mat,1);
            
        end
        
        %Compute EEG variance
        function varMat = getEEGVar(obj,subject,condition,stimulation, bLaplace)
            
            mat = obj.SelectMatrix(subject,condition,stimulation, bLaplace);
            varMat = var(mat,[],1);
            
        end
        
        %Normalize EEG
        function NEEG = NormalizeEEG(obj,subject,condition,stimulation,varargin)
            
            [results] = parseInput(obj.mode,varargin{:});
            
            %Select the matrix, EEG is also left in grouped cases
            [EEG,rEEG] = obj.SelectGroupMatrix(subject,condition,stimulation,1,results.Group);
            
            %Get baseline subtracted signal + variability
            for i = 1:size(EEG,1)
                mu = mean(EEG(i,1:251));
                sd = std(EEG(i,:));
                
                B(i,:) = EEG(i,:) - mu;
                S(:,i) = sd;
            end
            
            NEEG = mean(B,1)/mean(S);
            
            if (~isempty(rEEG))
                for i = 1:size(rEEG,1)
                    rmu = mean(rEEG(i,1:251));
                    rsd = std(rEEG(i,:));
                    
                    rB(i,:) = rEEG(i,:) - rmu;
                    rS(:,i) = rsd;
                end
                rNEEG = mean(rB,1)/mean(rS);
                NEEG = cat(3,NEEG,rNEEG);
            end
            
        end
        
        %Normalize EEGs across all subjects
        function EEGObj = NormalizeAllEEG(obj,condition,stimulation,varargin)
            
            EEGObj.data = [];
            for i = 1:size(obj.data,3)
                EEGObj.data = cat(1,EEGObj.data,obj.NormalizeEEG(i,condition,stimulation,varargin{:}));
            end
            
            EEGObj.condition = condition;
            EEGObj.stimulation = stimulation;
            
            if (size(EEGObj.data,3) == 2)
                EEGObj.lr = 1;
            else
                EEGObj.lr = 0;
            end
            
        end
        
        %Plot subjects individually
        function hFig = EEGAllZFig(obj,EEGObj)
            
            c_lr = {'Left','Right'};
            for lr = 1:size(EEGObj.data,3)
                if (EEGObj.lr)
                    d_type = c_lr{lr};
                else
                    d_type = 'Ungrouped';
                end
                for i = 1:size(EEGObj.data,1)
                    meanMat = EEGObj.data(i,:,lr);
                    figure;
                    plot(obj.Time,meanMat,'b', 'LineWidth', 1.5);
                    hold on;
                    plot([obj.Time(1) obj.Time(end)],[0 0], 'k--', 'LineWidth', 1); hold on;
                    y = [-0.1 0.1];
                    plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5);
                    hFig = gcf;
                    title(['Subject ' num2str(i) ' ' d_type]);
                end
            end
        end
        
        %Compute differences between left and right
        function EEGO = GetEEGDiff(~,EEGObj)
            EEGO = EEGObj;
            if (size(EEGObj.data,3) > 1)
                EEGO.data = EEGObj.data(:,:,1) - EEGObj.data(:,:,2);
            else
                disp('Object contains ungrouped data...');
            end
        end
                
        %Compute an average between 2 EEG data sets (Example, anodalpre
        %averaged with cathodalpre)
        function avgEEGObj = Avg2Sets(~,EEGObj1,EEGObj2)
            
            z = (EEGObj1.data + EEGObj2.data)./2;
            
            avgEEGObj.data = z;
            avgEEGObj.condition = 1;
            avgEEGObj.stimulation = 3;
            
        end
        
        %Compute difference between 2 EEG data sets (Example, post-pre)
        function diffEEGObj = Diff2Sets(~,EEGObj2,EEGObj1)
            
            z = EEGObj2.data - EEGObj1.data;
            
            diffEEGObj.data = z;
            diffEEGObj.condition = 3;
            diffEEGObj.stimulation = EEGObj2.stimulation;
        end
        
        %Significant from zero testing
        function pTest = TestEEG(~,EEGObj)
            
            for i = 1:size(EEGObj.data,3)
                for j = 1:size(EEGObj.data,2)
                    [~,pTest(i,j)] = ttest(EEGObj.data(:,j,i));
                end
            end
        end
        
        %Plotting Function for unsmoothed EEG from data 
        function hFig = EEGFig(obj,EEGObj,pTest)
            
            meanMat = mean(EEGObj.data,1);
            varMat = var(EEGObj.data,[],1);
            
            posVarMat = meanMat + varMat;
            negVarMat = meanMat - varMat;
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            figure;
            plot(obj.Time,meanMat,'k', 'LineWidth', 1.5);
            hold on;
            plot([obj.Time(1) obj.Time(end)],[0 0], 'k--', 'LineWidth', 1); hold on;
            patch([obj.Time fliplr(obj.Time)], [posVarMat fliplr(negVarMat)],...
                'k', 'FaceAlpha', 0.3, 'EdgeColor', 'w');
            
            y = [-0.5 0.5];
            plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5); hold on;
            xlim([X_MIN X_MAX]);
            ylim([y(1) y(2)]);
            if (obj.a_type == 0)
                %Include patch for 90% movement
                timeNP = obj.Time(obj.nP(1)+251:obj.nP(2)+251);
                patch([timeNP fliplr(timeNP)], [y(2)*ones(length(timeNP),1)' y(1)*ones(length(timeNP),1)'], 'c',...
                    'FaceAlpha',0.1,'EdgeColor','c');
            end
            hFig = gcf;
            
            
            %If pTest input given then add a significance patch
            if (nargin == 3)
                hSig = find(pTest < 0.05);
                ranges = group(hSig,1,3);
                
                %Now for each group
                for i = 1:length(ranges)/2
                    %Patch the area red
                    patchArea = hSig(ranges(2*i-1):ranges(2*i));
                    patch([obj.Time(patchArea) ...
                        fliplr(obj.Time(patchArea))], ...
                        [posVarMat(patchArea) fliplr(negVarMat(patchArea))], 'r', ...
                        'FaceAlpha', 0.5, 'EdgeColor', 'r');
                    
                end
            end
            
        end
        
        %Plotting function that automates the ptest and plotting
        function testedFig = PlotRawTest(obj,EEGObj)
            
            p = obj.TestEEG(EEGObj); 
            testedFig = obj.EEGFig(EEGObj,p); 
            
        end
        
        %Returns sliding window average EEG (24 sample default window)
        function sEEG = SWinEEG(~,EEGObj)
            winSize = 24;
            eegdata = EEGObj.data;
            
            %First extract variance of signal
            
            for i = winSize/2+1:size(eegdata,2)-winSize/2
                sAvg(:,i) = mean(eegdata(:,i-winSize/2:i+winSize/2),2);
            end
            
            EEGObj.data = sAvg;
            sEEG = EEGObj;
            sEEG.winSize = 24/2;
        end
        
        %Plot sliding windowed EEG with significance testing if given
        function hFig = PlotSWin(obj,EEGObj,pTest)
            
            meanMat = mean(EEGObj.data,1);
            varMat = var(EEGObj.data,[],1);
            
            posVarMat = meanMat + varMat;
            negVarMat = meanMat - varMat;
            winSize = EEGObj.winSize;
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            figure;
            plot(obj.Time(winSize+1:end-winSize),meanMat(winSize+1:end),'k', 'LineWidth', 1.5);
            hold on;
            plot([obj.Time(winSize+1) obj.Time(end-winSize)],[0 0], 'k--', 'LineWidth', 1); hold on;
            p = patch([obj.Time(winSize+1:end-winSize) fliplr(obj.Time(winSize+1:end-winSize))], [posVarMat(winSize+1:end) fliplr(negVarMat(winSize+1:end))],...
                'k', 'FaceAlpha', 0.3, 'EdgeColor', 'w');
            uistack(p,'bottom'); 
            %Set limits properly
            y = [-0.5 0.5];
            plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5); hold on;
            xlim([X_MIN X_MAX]);
            ylim([-0.5 0.5]);
            
            
            
            %If pTest input given then add a significance patch
            if (nargin == 3)
                hSig = find(pTest < 0.05);
                ranges = group(hSig,1,3);
                
                %Now for each group
                for i = 1:length(ranges)/2
                    %Patch the area red
                    patchArea = hSig(ranges(2*i-1):ranges(2*i));
                    patch([obj.Time(patchArea) ...
                        fliplr(obj.Time(patchArea))], ...
                        [posVarMat(patchArea) fliplr(negVarMat(patchArea))], 'r', ...
                        'FaceAlpha', 0.5, 'EdgeColor', 'r');
                    
                end
            end
            
            if (obj.a_type == 0)
                %Include patch for 90% movement
                timeNP = obj.Time(obj.nP(1)+251:obj.nP(2)+251);
                patch([timeNP fliplr(timeNP)], [y(2)*ones(length(timeNP),1)' y(1)*ones(length(timeNP),1)'], 'c',...
                    'FaceAlpha',0.1,'EdgeColor','c');
            end
            
            hFig = gcf;
            
        end
        
        %Function that processes smoothing and outputs a figure 
        function smoothPlot = PlotSmoothed(obj,EEGObj)
            sEEG = obj.SWinEEG(EEGObj);
            pEEG = obj.TestEEG(sEEG);
            smoothPlot = obj.PlotSWin(sEEG,pEEG);
            
        end
        
        %Compute all statistics for a selected EEG window
        function windowEEG = PickWin(~,EEGObj,winVector)
            windowEEG = EEGObj;
            windowEEG.data = EEGObj.data(:,winVector);
            windowEEG.winAvg = mean(windowEEG.data,2);
            [windowEEG.hTest, windowEEG.pTest] = ttest(windowEEG.winAvg);
            windowEEG.meanAvg = mean(windowEEG.winAvg);
            windowEEG.varAvg = var(windowEEG.winAvg);
            windowEEG.stdEr = std(windowEEG.winAvg)/sqrt(length(windowEEG.winAvg));
        end
               
        
        %% Behavioural Correlations
        
        %Grab a matrix for a given type of behavioural grouping
        function specMat = SpecifyMatrix(obj,subject,condition,stimulation,group)
            
            %Select conditions and get laplacian
            selMat = vertcat(obj.data{stimulation,condition,subject});
            selMat = selMat - vertcat(obj.channels{stimulation,condition,subject});
            
            %Specify the conditionals of the request
            subInd = ismember(obj.T.Subject,subject);
            stimInd = ismember(obj.T.Stim,stimulation);
            condInd = ismember(obj.T.Cond,condition);
            selInd = subInd & stimInd & condInd;
            
            %Define new table containing values for specific Sub/Cond/Stim
            selT = obj.T(selInd,:);
            
            %Check if RT is a group 
            if sum(strcmpi(group,'RT'))
                bins = 2; 
                %Temporary indexer for RT Grand Median Split
                selT.MEDIAN = (selT.RT > obj.medianTime) + 1; 
                %RT for within subject median split 
%                 selT.MEDIAN = (selT.RT > median(selT.RT)) + 1; 
                group(strcmpi(group,'RT')) = {'MEDIAN'}; 
            end
            str_statement = '';  
            for i = length(group):-1:1
                
                d_val{i}  = 1 : max(unique(selT.(group{i})));
                d_val{i}(d_val{i} == 0 | isnan(d_val{i})) = [];
                str_statement = [str_statement ',' 'd_val{' num2str(i) '}']; 
                col_ids(length(group) - i + 1) = find(strcmpi([obj.T.Properties.VariableNames 'MEDIAN'],group{i})); 
                
            end
            str_statement(1) = []; 
            sets = eval(['allcomb(' str_statement ')']);             
            
            selCols = selT{:,col_ids}; 
            
            
            %From table group the trials into a cell
            for s = 1 : size(sets,1) 
                specMat{s} = selMat(ismember(selCols,sets(s,:),'rows'),:);
                
                %If trial deficient, remove to prevent noisy weighting,
                %mainly for grand average median split 
                if (sum(ismember(selCols,sets(s,:),'rows')) < 10)
                    specMat{s} = []; 
                end
                
            end
        end
        
        %Compute normalized EEG blocks for each behavioural grouping
        function specEEG = ComputeSpecEEG(obj,subject,condition,stimulation,varargin)
            
            p = inputParser;
            
            defaultGroup = 'Target';
            validGroup = [obj.T.Properties.VariableNames(6:end),'Combine'];
            checkGroup = @(x) any(validatestring(x,validGroup));
            addParameter(p, 'Group', defaultGroup, checkGroup);
            p.KeepUnmatched = true;
            parse(p,varargin{1:2});
            
            %Grab data
            if sum(strcmpi(p.Results.Group,{'Combine','RT'})) == 0
                specMat = obj.SpecifyMatrix(subject,condition,stimulation,{p.Results.Group}); 
            elseif strcmpi(p.Results.Group,'Combine') 
                specMat = obj.SpecifyMatrix(subject,condition,stimulation,varargin{3:end}); 
            else
                if length(varargin) ~= 3
                    numBins = 5; 
                else
                    numBins = varargin{3}; 
                end
                specMat = obj.SelBinMat(subject,condition,stimulation,numBins,true);
            end
            
            %Normalize
            for k =1:length(specMat)
                
                EEG = specMat{k};
                
                B = []; 
                S = []; 
                for i = 1:size(EEG,1)
                    mu = mean(EEG(i,1:251));
                    sd = std(EEG(i,:));
                    
                    B(i,:) = EEG(i,:) - mu;
                    S(:,i) = sd;
                end
                
                specEEG{k} = mean(B,1)/mean(S);
            end
        end
        
        %Compute EEG blocks for each behavioural grouping across all
        %subjects
        function EEGArray = ComputeAllSpecEEG(obj,condition,stimulation,varargin)
            for i = 1:size(obj.data,3)
                EEGArray(i,:) = obj.ComputeSpecEEG(i,condition,stimulation,varargin{:});
            end
        end
        
        %Compute TFR, grouped behaviourally
        function TFRArray = ComputeSpecTFR(obj,subject,condition,stimulation,varargin)
            
            p = inputParser;
            
            defaultGroup = 'Target';
            validGroup = obj.T.Properties.VariableNames(6:end);
            checkGroup = @(x) any(validatestring(x,validGroup));
            addParameter(p, 'Group', defaultGroup, checkGroup);
            
            p.KeepUnmatched = true;
            parse(p,varargin{1:2});
            
            %Grab data
            if (~strcmpi(p.Results.Group,'RT'))
                specMat = obj.SpecifyMatrix(subject,condition,stimulation,{p.Results.Group});
            else
                %If RT return binned data (default bin = 5, laplacian =
                %true)
                
                %Prompt user for current instance's RT
                %If third varargin provided use that here, otherwise use 5
                if length(varargin) ~= 3
                    numBins = 5; 
                else
                    numBins = varargin{3}; 
                end
                specMat = obj.SelBinMat(subject,condition,stimulation,numBins,true);
            end
            
            %Conditionals (Using Pseudo-Z) 
            TFRMat = cell(1,length(specMat)); 
            for i = 1:length(specMat)
                TFRMat{i} = MakeTFRP(specMat{i}',obj.Freq,obj.Time,obj.Fs,7,[obj.Time(1) obj.Time(251)], 1);
            end
            for i = 1 : length(specMat)
                TFRArray(1,i) = TFRObject(TFRMat{i}, 1, subject, condition, stimulation, 1);
            end
        end
        
        %Compute TFR, grouped behaviourally, all subjects
        function TFRArray = ComputeAllSpecTFR(obj,condition,stimulation,varargin)
            
            TFRCell = cell(1,12);
            tic;
            parfor i = 1:size(obj.data,3)
                TFRCell{i} = obj.ComputeSpecTFR(i,condition,stimulation,varargin{:});
            end
            
            for i = 1 : length(TFRCell) 
                TFRArray(i,:) = TFRCell{i}; 
            end
            
            toc;
            
        end
        
        %Perform a block subtraction across two sets (TFR1 - TFR2)
        function subArray = SubtractBlock(~,TFR1,TFR2)
            
            for i = 1:numel(TFR1)
                
                subTFR = TFR1(i).GetTFR - TFR2(i).GetTFR;
                subArray(i) = TFR1(i);
                subArray(i).TFR = subTFR;
                
            end
            
            %Reform the original matrix
            subArray = subArray';
            subArray = reshape(subArray,size(TFR1));
            
        end
        
        %Regress within a given TFR across TFR values (depreciated, use
        %RegWinTFR for access to frequency windows)
        function regTFR = RegSpecTFR(obj,TFRBlock,regMat,win)
            
            if nargin < 4
                win = 1:length(obj.Time);
            end
            
            regMat = [ones(length(regMat),1) regMat'];
            
            for sub = 1:size(TFRBlock,1)
                
                %Extract target TFRs
                for targ = 1:size(TFRBlock,2)
                    
                    TFR(:,:,targ) = TFRBlock(sub,targ).GetTFR;
                    
                end
                
                %Perform regressions across third dimension with feature
                %input regMat
                for j = win(1):win(end)
                    for k = 1:size(TFR,1)
                        b = regress(squeeze(TFR(k,j,:)),regMat);
                        beta(k,j) = b(2);
                    end
                end
                
                %Create coeff matrix and meet sizing specification
                regTFR(:,:,sub) = [beta zeros(k,length(j+1:size(TFR,2)))];
                
            end
            
        end
        
        %Regress within a given TFR window across TFR values
        function regTFR = RegWinTFR(~,TFRBlock,regMat,win,fwin)
            
            regMat = [ones(length(regMat),1) regMat'];
            
            for sub = 1:size(TFRBlock,1)
                
                %Extract target TFRs
                for targ = 1:size(TFRBlock,2)
                    
                    TFR(:,:,targ) = TFRBlock(sub,targ).GetTFR;
                    
                end
                
                %Take a windowed average for variable
                avgTFR = squeeze(mean(mean(TFR(fwin,win,:))));
                
                %Perform regressions across third dimension with feature
                %input regMat
                for i = 1:length(avgTFR)
                    b = regress(avgTFR,regMat);
                    regTFR(sub) = b(2);
                end
                
                %Create coeff matrix and meet sizing specification
                
                
            end
            
        end
        
        %Regress within a given EEG window
        function regEEG = RegSpecEEG(~,EEGBlock,regMat,win)
            
            regMat = [ones(length(regMat),1) regMat'];
            
            %Compute average of each window in the EEGBlock
            for sub = 1:size(EEGBlock,1)
                
                subBlock = EEGBlock(sub,:);
                
                for i = 1:size(EEGBlock,2)
                    specBlock = subBlock{i};
                    avgEEG(i) = mean(specBlock(win));
                end
                
                %Now regress
                b = regress(avgEEG',regMat);
                regEEG(sub) = b(2);
            end
        end
        
        %Correlate with input data, slope TFR is given by regTFR (slopes
        %across TFRs across trial conditions)
        function [corrTFR,betaTFR] = CorrSpecTFR(obj,SlopeTFR,inputData,win)
            
            if nargin < 4
                win = 1:length(obj.Time);
            end
            
            if nargin < 5
                thres = 1;
            end
            
            %If it's a wrapped object extract TFR, otherwise if input is a
            %straight TFR just transfer over
            if(isa(SlopeTFR(1),'TFRObject'))
                for sub = 1:length(SlopeTFR)
                    TFR(:,:,sub) = SlopeTFR(sub).GetTFR;
                end
                
                SlopeTFR = TFR;
            end
            
            %Do a windowed analysis of the data given by win
            for k = win(1):win(end)
                for j = 1:size(SlopeTFR,1)
                    
                    [rscore, P] = corrcoef(squeeze(SlopeTFR(j,k,:))',inputData);
                    corrTFR(j,k) = rscore(2);
                    thresTFR(j,k) = P(2);
                    
                    %Get beta values
                    fitMat = polyfit(inputData(:),squeeze(SlopeTFR(j,k,:)),1);
                    betaTFR(j,k) = fitMat(1);
                end
            end
            
            corrTFR = [corrTFR zeros(j,length(k+1:size(SlopeTFR,2)))];
            betaTFR = [betaTFR zeros(j,length(k+1:size(SlopeTFR,2)))];
            thresTFR = [thresTFR nan(j,length(k+1:size(SlopeTFR,2)))];
            
            %Clean up thresholding matrix
            thresTFR(:,1:win(1)-1) = NaN;
            
            if (thres)
                %Apply threshold
                corrTFR(thresTFR > 0.05) = NaN;
                betaTFR(thresTFR > 0.05) = NaN;
            end
        end
        
        %Extract window components of TFR
        function waTFR = ExtractTFRAvg(obj,TFRArray,winTime,winFreq)
            
            %Get time window in samples
            rawSpecWindowT = 251 + obj.baseOffset + round((winTime).*obj.Fs./1000); 
            
            %Fill in gaps for time window
            specWindowT = rawSpecWindowT(1):rawSpecWindowT(end); 
            specWindowF = winFreq; 
            
            %Mean selection function
            GetWindowMean = @(TFR) (mean(mean(TFR.TFR(specWindowF,specWindowT)))); 
            
            %Mean computation 
            waTFR = arrayfun(GetWindowMean,TFRArray); 
                        
        end
        
        %Extract window components of EEG
        function waEEG = ExtractEEGAvg(obj,specEEG,winTime)
            
            %Get the sample window of interest
            rawSpecWindow = 251 + obj.baseOffset + round((winTime).*obj.Fs./1000);
            
            %Fill in the gaps
            specWindow = rawSpecWindow(1):rawSpecWindow(end);
            
            %For each cell, extract out the window and compute the average
            waEEG = cellfun(@mean,cellfun(@(x) x(specWindow), specEEG,'UniformOutput',false));
            
            
        end
        
        %Quick plot function for extracted windowed averages
        function wFig = PlotGroupedWinEEG(obj,specEEG,winTime)
            
            waEEG = obj.ExtractEEGAvg(specEEG,winTime);
            figure;
            plot(waEEG','.', 'MarkerSize', 12,'Color', 'k'); hold on;
            plot(mean(waEEG,1),'s','MarkerSize',5,'Color','r','MarkerFaceColor','r');
            xlabel('Target');
            ylabel('Normalized Average EEG');
            title([num2str(winTime(1)) '-' num2str(winTime(end)) 'ms']);
            set(gca, 'XTick', 1:size(waEEG,2));
            xlim([0 size(waEEG,2)+1]);
            ylim([-1 1]);
            wFig = gcf;
            
        end
        
        %Perform comparison dependent on group type (EEG)
        function regEEG = RegGroupedEEG(obj,specEEG,winTime,varargin)
            
            %Do some parsing to select the grouping type
            results = parseInput(obj.mode,varargin{:});
            
            %Compute the windowed averages for each subject and group
            %bin
            waEEG = obj.ExtractEEGAvg(specEEG,winTime);
            
            %Once the group type is selected perform the appropriate
            %regressions/subtractions
            regEEG = ExtractRegMat(waEEG,results.Group);
            
        end
        
        %Perform comparison dependent on group type (TFR)
        function regTFR = RegGroupedTFR(obj,specTFR,winTime,winFreq,varargin) 
            
            results = parseInput(obj.mode,varargin{:}); 
            waTFR = obj.ExtractTFRAvg(specTFR,winTime,winFreq); 
            regTFR = ExtractRegMat(waTFR,results.Group); 
            
        end
        
        %Significant from zero testing
        function pTest = CTestEEG(~,cell_EEG)
            
            %Collapse each cell column vector into array stacks
            for i = 1 : size(cell_EEG,2)
                stacked{i} = vertcat(cell_EEG{:,1});
            end
            %Perform a ttest for each of the cells 
            pTest = cellfun(@(x) ttest(x),stacked,'un',0);       
        end
        
        %Plotting Function for unsmoothed EEG from data 
        function hFig = CEEGFig(obj,cell_EEG,pTest)
            
            %Stack the EEG
            EEG = vertcat(cell_EEG{:}); 
            
            meanMat = mean(EEG,1);
            varMat = var(EEG,[],1);
            
            posVarMat = meanMat + varMat;
            negVarMat = meanMat - varMat;
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            figure;
            plot(obj.Time,meanMat,'k', 'LineWidth', 1.5);
            hold on;
            plot([obj.Time(1) obj.Time(end)],[0 0], 'k--', 'LineWidth', 1); hold on;
            p = patch([obj.Time fliplr(obj.Time)], [posVarMat fliplr(negVarMat)],...
                'k', 'FaceAlpha', 0.3, 'EdgeColor', 'w');
            uistack(p,'bottom'); 
            
            y = [-0.5 0.5];
            plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5); hold on;
            xlim([X_MIN X_MAX]);
            ylim([y(1) y(2)]);
            if (obj.a_type == 0)
                %Include patch for 90% movement
                timeNP = obj.Time(obj.nP(1)+251:obj.nP(2)+251);
                patch([timeNP fliplr(timeNP)], [y(2)*ones(length(timeNP),1)' y(1)*ones(length(timeNP),1)'], 'c',...
                    'FaceAlpha',0.1,'EdgeColor','c');
            end
            hFig = gcf;
            
            
            %If pTest input given then add a significance patch
            if (nargin == 3)
                hSig = find(pTest);
                ranges = group(hSig,1,3);
                
                %Now for each group
                for i = 1:length(ranges)/2
                    %Patch the area red
                    patchArea = hSig(ranges(2*i-1):ranges(2*i));
                    patch([obj.Time(patchArea) ...
                        fliplr(obj.Time(patchArea))], ...
                        [posVarMat(patchArea) fliplr(negVarMat(patchArea))], 'r', ...
                        'FaceAlpha', 0.5, 'EdgeColor', 'r');
                    
                end
            end
            
        end
        
        %Plotting function that automates the ptest and plotting
        function testedFig = CPlotRawTest(obj,cell_EEG)
            
            p = obj.CTestEEG(cell_EEG); 
            testedFig = obj.CEEGFig(cell_EEG,p{1}); 
            
        end
        
        %Returns sliding window average EEG (24 sample default window)
        function sEEG = CSWinEEG(~,cell_EEG)
            winSize = 24;
            eegdata = vertcat(cell_EEG{:});
            
            %First extract variance of signal
            
            for i = winSize/2+1:size(eegdata,2)-winSize/2
                sAvg(:,i) = mean(eegdata(:,i-winSize/2:i+winSize/2),2);
            end
            
            EEGObj.data = sAvg;
            sEEG = EEGObj;
            sEEG.winSize = 24/2;
        end
        
        %Plot sliding windowed EEG with significance testing if given
        function hFig = CPlotSWin(obj,cell_EEG,pTest)
            
            EEG = vertcat(cell_EEG{:}); 
            
            meanMat = mean(EEG,1);
            varMat = var(EEG,[],1);
            
            posVarMat = meanMat + varMat;
            negVarMat = meanMat - varMat;
            winSize = EEGObj.winSize;
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            
            figure;
            plot(obj.Time(winSize+1:end-winSize),meanMat(winSize+1:end),'k', 'LineWidth', 1.5);
            hold on;
            plot([obj.Time(winSize+1) obj.Time(end-winSize)],[0 0], 'k--', 'LineWidth', 1); hold on;
            patch([obj.Time(winSize+1:end-winSize) fliplr(obj.Time(winSize+1:end-winSize))], [posVarMat(winSize+1:end) fliplr(negVarMat(winSize+1:end))],...
                'k', 'FaceAlpha', 0.3, 'EdgeColor', 'w');
            
            %Set limits properly
            y = [-0.5 0.5];
            plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5); hold on;
            xlim([X_MIN X_MAX]);
            ylim([-0.5 0.5]);
            
            
            
            %If pTest input given then add a significance patch
            if (nargin == 3)
                hSig = find(pTest < 0.05);
                ranges = group(hSig,1,3);
                
                %Now for each group
                for i = 1:length(ranges)/2
                    %Patch the area red
                    patchArea = hSig(ranges(2*i-1):ranges(2*i));
                    patch([obj.Time(patchArea) ...
                        fliplr(obj.Time(patchArea))], ...
                        [posVarMat(patchArea) fliplr(negVarMat(patchArea))], 'r', ...
                        'FaceAlpha', 0.5, 'EdgeColor', 'r');
                    
                end
            end
            
            if (obj.a_type == 0)
                %Include patch for 90% movement
                timeNP = obj.Time(obj.nP(1)+251:obj.nP(2)+251);
                p = patch([timeNP fliplr(timeNP)], [y(2)*ones(length(timeNP),1)' y(1)*ones(length(timeNP),1)'], 'c',...
                    'FaceAlpha',0.1,'EdgeColor','c');
            end
            
            hFig = gcf;
            
        end
        
        %Function that processes smoothing and outputs a figure 
        function smoothPlot = CPlotSmoothed(obj,cell_EEG)
            sEEG = obj.CSWinEEG(cell_EEG);
            pEEG = obj.TestEEG(sEEG);
            smoothPlot = obj.PlotSWin(sEEG,pEEG);
            
        end
        
        %Compute all statistics for a selected EEG window
        function windowEEG = CPickWin(~,cell_EEG,winVector)
            
            EEG = vertcat(cell_EEG{:}); 
            
            windowEEG = EEG;
            windowEEG.data = EEG(:,winVector);
            windowEEG.winAvg = mean(windowEEG.data,2);
            [windowEEG.hTest, windowEEG.pTest] = ttest(windowEEG.winAvg);
            windowEEG.meanAvg = mean(windowEEG.winAvg);
            windowEEG.varAvg = var(windowEEG.winAvg);
            windowEEG.stdEr = std(windowEEG.winAvg)/sqrt(length(windowEEG.winAvg));
        end
        
        %Generate an overlayed version of comparisons time-series 
        function hFig = OverlayEEGStats(obj,EEG,smoothOpt,legname)
            
            %Filter out empties 
            EEG1Empt = cellfun(@isempty,EEG(:,1),'un',0); 
            EEG2Empt = cellfun(@isempty,EEG(:,2),'un',0);
            EEGEmpt = sum(cell2mat([EEG1Empt,EEG2Empt]),2) > 0; 
            EEG(EEGEmpt,:) = []; 
            
            %Take difference
            D = cellfun(@(x,y) x - y, EEG(:,1), EEG(:,2),'un',0);
            c1 = hex2rgb('ef8a62')./255 - 10/255; 
            c2 = hex2rgb('67a9cf')./255 - 10/255; 
            
            %Define CPP window (Kelly and Connell, 2013) 
%             cpp_win = [200,350;-250,-100]; 
            
            figure;
            switch smoothOpt
                
                case 0
                    p = obj.CTestEEG(D); %Compute stats on difference
                    pTest = p{1}; 
                    p1 = obj.PlotOverlayEEGStats(EEG(:,1)); 
                    set(p1(1),'Color',c1);
                    set(p1(2),'FaceColor',c1); 
                    p2 = obj.PlotOverlayEEGStats(EEG(:,2));
                    set(p2(1),'Color',c2); 
                    set(p2(2),'FaceColor',c2);
                case 1
                    sEEG = obj.CSWinEEG(D); 
                    p = obj.TestEEG(sEEG); 
                    pTest = p < 0.05; 
                    sEEG1 = obj.CSWinEEG(EEG(:,1)); 
                    sEEG2 = obj.CSWinEEG(EEG(:,2)); 
                    p1 = obj.PlotSmoothOverlayEEGStats(sEEG1); 
                    set(p1(1),'Color',c1); 
                    set(p1(2),'FaceColor',c1); 
                    p2 = obj.PlotSmoothOverlayEEGStats(sEEG2); 
                    set(p2(1),'Color',c2); 
                    set(p2(2),'FaceColor',c2); 
                    
                    
            end
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            y = [-0.5 0.5];
            plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5); hold on;
            xlim([X_MIN X_MAX]);
            ylim([y(1) y(2)]);
            
            plot([obj.Time(1) obj.Time(end)],[0 0], 'k--', 'LineWidth', 1); hold on;
            if (obj.a_type == 0)
                %Include patch for 90% movement
                timeNP = obj.Time(obj.nP(1)+251:obj.nP(2)+251);
               pm = patch([timeNP fliplr(timeNP)], [y(2)*ones(length(timeNP),1)' y(1)*ones(length(timeNP),1)'], 'c',...
                    'FaceAlpha',0.1,'EdgeColor','none');
                uistack(pm,'bottom');
            end
            
            %Significance patches
            hSig = find(pTest);
            ranges = group(hSig,1,3);
            
            %Now for each group
            for i = 1:length(ranges)/2
                %Patch the area red
                patchArea = hSig(ranges(2*i-1):ranges(2*i));
                pp = patch([obj.Time(patchArea) ...
                    fliplr(obj.Time(patchArea))], ...
                    [ones(size(patchArea)).*y(2) ones(size(patchArea)).*y(1)], [0.6 0.6 0.6], ...
                    'FaceAlpha', 0.5, 'EdgeColor', 'none');
                uistack(pp,'bottom'); 
            end 
            
            %CPP Window 
%             plotvline(cpp_win(obj.a_type+1,:)); 
            
            legend([p1(1) p2(1)],[legname(1),legname(2)]); 
            
        end
        
        function hFig = OverlayEEG(obj,EEG,smoothOpt,legname)
            
            %Take difference
            D = cellfun(@(x,y) x - y, EEG(:,1), EEG(:,2),'un',0);
            c1 = hex2rgb('ef8a62')./255 - 10/255; 
            c2 = hex2rgb('67a9cf')./255 - 10/255; 
            figure;
            
            switch smoothOpt
                case 0
                    %Create an overlayed plot
                    p1 = obj.PlotOverlayEEGStats(EEG(:,1)); 
                    set(p1(1),'Color',c1);
                    set(p1(2),'FaceColor',c1); 
                    p2 = obj.PlotOverlayEEGStats(EEG(:,2));
                    set(p2(1),'Color',c2); 
                    set(p2(2),'FaceColor',c2);
                case 1
                    sEEG1 = obj.CSWinEEG(EEG(:,1)); 
                    sEEG2 = obj.CSWinEEG(EEG(:,2)); 
                    p1 = obj.PlotSmoothOverlayEEGStats(sEEG1); 
                    set(p1(1),'Color',c1); 
                    set(p1(2),'FaceColor',c1); 
                    p2 = obj.PlotSmoothOverlayEEGStats(sEEG2); 
                    set(p2(1),'Color',c2); 
                    set(p2(2),'FaceColor',c2); 
            end
            
            xl = [-0.5 -0.8];
            X_MIN = xl(obj.a_type + 1);
            X_MAX = 1.0;
            y = [-0.5 0.5];
            plot([0 0], [y(1) y(2)], 'k:', 'LineWidth', 1.5); hold on;
            xlim([X_MIN X_MAX]);
            ylim([y(1) y(2)]);
            
            plot([obj.Time(1) obj.Time(end)],[0 0], 'k--', 'LineWidth', 1); hold on;
            if (obj.a_type == 0)
                %Include patch for 90% movement
                timeNP = obj.Time(obj.nP(1)+251:obj.nP(2)+251);
               pm = patch([timeNP fliplr(timeNP)], [y(2)*ones(length(timeNP),1)' y(1)*ones(length(timeNP),1)'], 'c',...
                    'FaceAlpha',0.1,'EdgeColor','none');
                uistack(pm,'bottom');
            end
            
            legend([p1(1) p2(1)],[legname(1),legname(2)]); 
        end
            
        function p = PlotOverlayEEGStats(obj,cell_EEG) 
            
            EEG = vertcat(cell_EEG{:}); 
            meanMat = mean(EEG,1);
            varMat = var(EEG,[],1);
            posVarMat = meanMat + varMat;
            negVarMat = meanMat - varMat;
            p(1) = plot(obj.Time,meanMat,'k', 'LineWidth', 1.5);
            hold on;
            p(2) = patch([obj.Time fliplr(obj.Time)], [posVarMat fliplr(negVarMat)],...
                'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
           
        end
        
        function p = PlotSmoothOverlayEEGStats(obj,sEEG)

            winSize = sEEG.winSize; 
            meanMat = mean(sEEG.data,1);
            varMat = var(sEEG.data,[],1);
            posVarMat = meanMat + varMat;
            negVarMat = meanMat - varMat;
            p(1) = plot(obj.Time(winSize+1:end-winSize),meanMat(winSize+1:end),'k', 'LineWidth', 1.5);
            hold on;
            p(2) = patch([obj.Time(winSize+1:end-winSize) fliplr(obj.Time(winSize+1:end-winSize))], [posVarMat(winSize+1:end) fliplr(negVarMat(winSize+1:end))],...
                'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            
        end
        
    end
    
end

%% Local methods
%Return a specified field of struct
function value = getFieldi(S,field)
names = fieldnames(S);
isField = strcmpi(field,names);

if any(isField)
    value = S.(names{isField});
else
    value = [];
end
end

%Input parsing function
function [results] = parseInput(eMode,varargin)
p = inputParser;

defaultType = 'Pseudo-Z';
validTypes = {'Pseudo-Z', 'Power'};
checkType = @(x) any(validatestring(x,validTypes));

defaultLaplace = 'Laplacian';
validLaplace = {'None', 'Laplacian'};
checkLaplace = @(x) any(validatestring(x,validLaplace));

defaultGroup = 'None';
if (eMode == 1) %RDM
    validGroup = {'L','R','LR','CW','RT'}; 
elseif (eMode == 2) %ARMREACH
    validGroup = {'IHP' 'None' 'Target' 'Reach' 'RHemifield' 'LHemifield',...
        'EHemifield' 'MedianSplit','TIHP','TargIHP','ReachVec'};
elseif (eMode == 3) %FC
    validGroup = {};
end
checkGroup = @(x) any(validatestring(x,validGroup));

addParameter(p, 'TFRtype', defaultType, checkType);
addParameter(p, 'SpatialFilter', defaultLaplace, checkLaplace);
addParameter(p, 'Group', defaultGroup, checkGroup);

p.KeepUnmatched = true;

parse(p,varargin{:});

gVal = find(strcmpi(validGroup,p.Results.Group));
results = p.Results;
end

function regR = ExtractRegMat(waEEG, results)

c_group = {'TIHP','targIHP','reachVec','EHemifield', ...
    'LHemifield','RHemifield','Reach','Target','IHP'};

%Check which grouping is being used and extract matching number
groupInd = find(strcmpi(c_group,results));

%Run conditionals and return compressed regression matrix for behaviours
if (groupInd == 2) %Special case for regression subtraction
    
    %Set up regressor
    cMat = [-10 -5 5 10];
    regMat = [ones(length(cMat),1) cMat'];
    
    for sub = 1 : size(waEEG,1)
        
        %Perform both regressions after converting waEEG to form that can be
        %regressed
        lReg = regress(waEEG(sub,1:4)',regMat);
        rReg = regress(waEEG(sub,5:end)',regMat);
        
        %Subtract regressions to give final output matrix for this condition
        regR(sub,1) = lReg(2) - rReg(2);
        
    end
    
elseif (groupInd == 3) %Special case for regression
    
    
    cMat = [-17.5 -12.5 -2.5 2.5 12.5 17.5]; 
    regMat = [ones(length(cMat),1) cMat']; 
    
    for sub = 1 : size(waEEG,1)
        
        %Perform target regression
        b = regress(waEEG(sub,:)',regMat);
        regR(sub,1) = b(2); 
    end
    
else %General case for subtraction (L-R) - combine all data here...
    %Average across all trials in this type of conditional... 
    
    regR = waEEG(:,1) - waEEG(:,2);
    
end

end



