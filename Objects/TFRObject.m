classdef TFRObject
    %TFROBJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c_range;
        type;
        condition;
        subject;
        TFR;
        stimulation;
        bLaplace;
        direction; 
    end
    
    methods
        
        function obj = TFRObject(TFR, type, subject, condition, stimulation, bLaplace,direction)
            c_types = [1 1; 0.1 0.0002]; 
            obj.type = type;
            obj.subject = subject;
            obj.condition = condition;
            obj.stimulation = stimulation;
            obj.c_range = c_types(type,bLaplace+1);
            obj.TFR = TFR;
            obj.bLaplace = bLaplace;
            
            if nargin < 7
                direction = []; 
            end
            
            obj.direction = direction; 
        end
        
        function [type] = GetType(obj)
            type = obj.type;
        end
        
        function subject = GetSubject(obj)
            subject = obj.subject;
        end
        
        function condition = GetCondition(obj)
            condition = obj.condition;
        end
        
        function c_range = GetRange(obj)
            c_range = obj.c_range;
        end
        
        function TFR = GetTFR(obj)
            TFR = obj.TFR;
        end
        
        function stimulation = GetStimulation(obj)
            stimulation = obj.stimulation;
        end
        
        function bLaplace = GetFilter(obj)
            bLaplace = obj.bLaplace; 
        end
        
        function direction = GetDirection(obj) 
            direction = obj.direction;
        end
    end
    
end

