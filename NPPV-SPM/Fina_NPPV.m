% ============================================================
% COPERNICUS MARINE SYSTEM (CMEMS)
% FINAL COMPLETE SCRIPT
% Cyclone Fina - NPPV
%
% Includes:
% 1. Mean map
% 2. Standard deviation map
% 3. Timeseries
% 4. Peak cyclone map
% 5. Linear regression trend
% 6. Anomaly vs Composite Baseline
% 7. Z-score anomaly vs Baseline
%
% IMPORTANT:
% Run BASELINE_NPPV script first and keep:
% baseline_mean
% baseline_std
%
% ============================================================

clear
clc
close all

%% ============================================================
%% TOOLBOXES
%% ============================================================

addpath 'C:\Programing\Digital Ocean\Report7abril\Exercise_Week6\Exercise_Week6\FinalProjectAlfred\m_map'
addpath 'C:\Programing\Digital Ocean\Report7abril\Exercise_Week6\Exercise_Week6\FinalProjectAlfred\unixtime2mat'

%% ============================================================
%% FOLDER
%% ============================================================

cd 'C:\Programing\Digital Ocean\Report7abril\Exercise_Week6\Exercise_Week6\FinalProjectAlfred'

output_dir = 'plots';

if ~exist(output_dir,'dir')
    mkdir(output_dir)
end

%% ============================================================
%% FILE
%% ============================================================

filename = 'NPP_Fina.nc';

ncdisp(filename)

%% ============================================================
%% READ VARIABLES
%% ============================================================

lon   = double(ncread(filename,'longitude'));
lat   = double(ncread(filename,'latitude'));
depth = double(ncread(filename,'depth'));

time = double(ncread(filename,'time'));
time = unixtime2mat(time);
time = datetime(time,'ConvertFrom','datenum');

npp = double(ncread(filename,'nppv'));

disp('--- ORIGINAL SIZE ---')
size(npp)

%% ============================================================
%% SURFACE LAYER
%% ============================================================

npp_surface = squeeze(npp(:,:,1,:));

%% ============================================================
%% MEAN / STD
%% ============================================================

npp_mean = mean(npp_surface,3,'omitnan');
npp_std  = std(npp_surface,0,3,'omitnan');

npp_mean = reshape(npp_mean,length(lon),length(lat));
npp_std  = reshape(npp_std ,length(lon),length(lat));

[Lon,Lat] = meshgrid(lon,lat);

%% ============================================================
%% MEAN MAP
%% ============================================================

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,npp_mean')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(jet)

title('Mean Surface NPPV - Cyclone Fina')

saveas(gcf,fullfile(output_dir,'NPP_Fina_mean.png'))

%% ============================================================
%% STD MAP
%% ============================================================

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,npp_std')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(parula)

title('STD Surface NPPV - Cyclone Fina')

saveas(gcf,fullfile(output_dir,'NPP_Fina_std.png'))

%% ============================================================
%% TIMESERIES
%% ============================================================

ts_npp = squeeze(mean(mean(npp_surface,1,'omitnan'),2,'omitnan'));

figure;

plot(time,ts_npp,'LineWidth',2)
grid on
box on

xlabel('Date')
ylabel('NPPV')
title('Cyclone Fina Daily Mean NPPV')

% Adaptive Y axis
ymax = max(ts_npp,[],'omitnan');
ymin = min(ts_npp,[],'omitnan');

padding = 0.10*(ymax-ymin);

if padding == 0
    padding = ymax*0.1 + 0.01;
end

ylim([max(0,ymin-padding) ymax+padding])

saveas(gcf,fullfile(output_dir,'NPP_Fina_timeseries.png'))

%% ============================================================
%% PEAK MAP
%% ============================================================

[~,idx_peak] = max(ts_npp,[],'omitnan');

peak_date = time(idx_peak);
peak_map  = npp_surface(:,:,idx_peak);

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,peak_map')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(jet)

title(['Peak NPPV Fina - ' datestr(peak_date)])

saveas(gcf,fullfile(output_dir,'NPP_Fina_peakmap.png'))

%% ============================================================
%% LINEAR REGRESSION
%% ============================================================

x = (1:length(ts_npp))';
y = ts_npp(:);

valid = ~isnan(y);
x = x(valid);
y = y(valid);

p = polyfit(x,y,1);
yfit = polyval(p,x);

SSres = sum((y-yfit).^2);
SStot = sum((y-mean(y)).^2);
R2 = 1 - SSres/SStot;

figure;

plot(x,y,'bo')
hold on
plot(x,yfit,'r','LineWidth',2)
grid on

xlabel('Time step')
ylabel('Mean NPPV')
title('Cyclone Fina Linear Trend')

legend('Observed','Linear Fit')

txt = sprintf('Slope = %.4f   R^2 = %.3f',p(1),R2);
text(min(x)+2,max(y)*0.95,txt)

saveas(gcf,fullfile(output_dir,'NPP_Fina_regression.png'))

%% ============================================================
%% ANOMALY VS BASELINE
%% ============================================================

if exist('baseline_mean','var')

fina_anomaly = npp_mean - baseline_mean;

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,fina_anomaly')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(jet)

title('Fina NPPV Anomaly vs Baseline')

saveas(gcf,fullfile(output_dir,'NPP_Fina_anomaly.png'))

end

%% ============================================================
%% Z SCORE VS BASELINE
%% ============================================================

if exist('baseline_mean','var') && exist('baseline_std','var')

fina_z = (npp_mean - baseline_mean) ./ baseline_std;

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,fina_z')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(parula)

title('Fina Standardized Anomaly (Z-score)')

saveas(gcf,fullfile(output_dir,'NPP_Fina_zscore.png'))

end