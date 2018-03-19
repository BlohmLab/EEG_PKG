function plotTFR( x, y, z, str, lab, xmin, xmax, t1, t2, Freq, c_range )
%PLOTTFR Summary of this function goes here
%   Detailed explanation goes here

figure('Name', str);
surf(x, y, z);
shading interp
view([0 0 1])
%plot3([0 0],[min(Freq) max(Freq)], 1.01*max(max(dTFR))*[1 1],'w:')
%caxis(max(abs(get(gca,'clim'))).*[-1 1])
caxis([-1*c_range c_range])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title([str 'Aligned to ' lab])
cc = colorbar;
ylabel(cc,'pseudo-Z score (re. baseline)')
colormap('Jet');
hold on;

xlim([xmin 1.0])

plot3([t1 t1], [0 max(Freq)], [1 1], 'k' ); 
plot3([t2 t2], [0 max(Freq)], [1 1], 'k' );

end

