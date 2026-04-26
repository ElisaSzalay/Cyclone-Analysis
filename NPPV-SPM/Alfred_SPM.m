% ============================================================
% COPERNICUS MARINE SYSTEM (CMEMS)
% FINAL COMPLETE SCRIPT
% Cyclone Alfred - SPM
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
% Run BASELINE_SPM script first and keep:
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

filename = 'SPM_Alfred.nc';

ncdisp(filename)

%% ============================================================
%% READ VARIABLES
%% ============================================================

lon = double(ncread(filename,'longitude'));
lat = double(ncread(filename,'latitude'));

time = double(ncread(filename,'time'));
time = unixtime2mat(time);
time = datetime(time,'ConvertFrom','datenum');

spm = double(ncread(filename,'SPM'));

disp('--- ORIGINAL SIZE ---')
size(spm)

% Expected:
% lon x lat x time

[Lon,Lat] = meshgrid(lon,lat);

%% ============================================================
%% MEAN / STD
%% ============================================================

spm_mean = mean(spm,3,'omitnan');
spm_std  = std(spm,0,3,'omitnan');

spm_mean = reshape(spm_mean,length(lon),length(lat));
spm_std  = reshape(spm_std ,length(lon),length(lat));

%% ============================================================
%% MEAN MAP
%% ============================================================

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,spm_mean')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(jet)

title('Mean Surface SPM - Cyclone Alfred')

saveas(gcf,fullfile(output_dir,'SPM_Alfred_mean.png'))

%% ============================================================
%% STD MAP
%% ============================================================

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,spm_std')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(parula)

title('STD Surface SPM - Cyclone Alfred')

saveas(gcf,fullfile(output_dir,'SPM_Alfred_std.png'))

%% ============================================================
%% TIMESERIES
%% ============================================================

ts_spm = squeeze(mean(mean(spm,1,'omitnan'),2,'omitnan'));

figure;

plot(time,ts_spm,'LineWidth',2)
grid on
box on

xlabel('Date')
ylabel('SPM (g m^{-3})')
title('Cyclone Alfred Daily Mean SPM')

% AUTO Y-AXIS ADAPTED TO DATA
ymax = max(ts_spm,[],'omitnan');
ymin = min(ts_spm,[],'omitnan');

padding = 0.10*(ymax-ymin);

if padding == 0
    padding = ymax*0.1 + 0.01;
end

ylim([max(0,ymin-padding) ymax+padding])

saveas(gcf,fullfile(output_dir,'SPM_Alfred_timeseries.png'))

%% ============================================================
%% PEAK MAP
%% ============================================================

[~,idx_peak] = max(ts_spm,[],'omitnan');

peak_date = time(idx_peak);
peak_map  = spm(:,:,idx_peak);

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

title(['Peak SPM Alfred - ' datestr(peak_date)])

saveas(gcf,fullfile(output_dir,'SPM_Alfred_peakmap.png'))

%% ============================================================
%% LINEAR REGRESSION
%% ============================================================

x = (1:length(ts_spm))';
y = ts_spm(:);

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
ylabel('Mean SPM')
title('Cyclone Alfred Linear Trend')

legend('Observed','Linear Fit')

txt = sprintf('Slope = %.4f   R^2 = %.3f',p(1),R2);
text(min(x)+2,max(y)*0.95,txt)

saveas(gcf,fullfile(output_dir,'SPM_Alfred_regression.png'))

%% ============================================================
%% ANOMALY VS BASELINE
%% ============================================================
% Requires baseline_mean loaded in workspace

if exist('baseline_mean','var')

alfred_anomaly = spm_mean - baseline_mean;

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,alfred_anomaly')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(jet)

title('Alfred SPM Anomaly vs Baseline')

saveas(gcf,fullfile(output_dir,'SPM_Alfred_anomaly.png'))

end

%% ============================================================
%% Z SCORE VS BASELINE
%% ============================================================

if exist('baseline_mean','var') && exist('baseline_std','var')

alfred_z = (spm_mean - baseline_mean) ./ baseline_std;

figure;

m_proj('mercator',...
'lon',[min(lon) max(lon)],...
'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,alfred_z')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12)

colorbar
colormap(parula)

title('Alfred Standardized SPM Anomaly (Z-score)')

saveas(gcf,fullfile(output_dir,'SPM_Alfred_zscore.png'))

end
