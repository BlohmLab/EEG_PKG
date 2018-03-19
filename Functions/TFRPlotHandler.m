classdef TFRPlotHandler
    %TFR plot handler. Defines TFR plot specifications...
    
    properties
        surfX;
        surfY;
        xMin;
        xMax
        Freq;
        tBaseStart;
        tBaseEnd;
    end
    
    methods
        function obj = TFRPlotHandler(surfx, surfy, xmin, xmax, freq, tbasestart, tbaseend)
            obj.surfX = surfx; 
            obj.surfY = surfy; 
            obj.xMin = xmin; 
            obj.xMax = xmax;
            obj.Freq = freq; 
            obj.tBaseStart = tbasestart; 
            obj.tBaseEnd = tbaseend; 
        end
        function plotTFR (obj, TFR, titleStr, c_range, optYlabel)
            %Plot TFR plots with saved instance fo plot handler with
            %set specifications
            if(nargin < 5 || isempty(optYlabel))
                optYlabel = 'pseudo-Z score (re. baseline)'; 
            end
           
            figure;
            surf(obj.surfX,obj.surfY, TFR);
            shading interp
            view([0 0 1])
            
            caxis([-1*c_range c_range])
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
            title(titleStr)
            cc = colorbar;
            ylabel(cc,optYlabel)
            colormap('Jet');
            hold on;
            
            xlim([obj.xMin obj.xMax])
            
            plot3([obj.tBaseStart obj.tBaseStart], [0 max(obj.Freq)], [1 1], 'k' );
            plot3([obj.tBaseEnd obj.tBaseEnd], [0 max(obj.Freq)], [1 1], 'k' );
             
            
            
        end
    end
    
end

