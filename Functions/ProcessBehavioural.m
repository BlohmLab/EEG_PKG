function [ selMeasure ] = ProcessBehavioural(datatable,meas_x,pol,grp)
%PROCESSBEHAVIOURAL

%Polarities 
polind = datatable.POLARITY == pol; 

%Group 
dirind = strcmpi(datatable.GROUP,grp); 

%Selection
selind = polind & dirind; 

%Return measure 
selMeasure = datatable.(upper(meas_x))(selind); 
end

