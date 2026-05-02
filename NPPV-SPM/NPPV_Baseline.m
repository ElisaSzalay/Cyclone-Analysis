% ============================================================
% COPERNICUS MARINE SYSTEM (CMEMS)
% FINAL BASELINE SCRIPT - NPPV COMPOSITE (3 PERIODS)
%
% Combines:
% baseline_NPPV_1.nc   (Nov 2019 - Mar 2020)
% baseline_NPPV_2.nc   (Nov 2020 - Mar 2021)
% baseline_NPPV_3.nc   (Nov 2021 - Mar 2022)
%
% Removes cyclone days
%
% Outputs:
% 1. Mean map
% 2. Standard deviation map
% 3. Daily timeseries
% 4. Peak event map
% 5. Linear regression trend
%
% Use this baseline to compare Alfred and Fina
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
%% FILE LIST
%% ============================================================

files = { ...
    'baseline_NPPV_1.nc',...
    'baseline_NPPV_2.nc',...
    'baseline_NPPV_3.nc'};

%% ============================================================
%% INITIALIZE
%% ============================================================

baseline_all  = [];
time_all      = [];

%% ============================================================
%% LOOP FILES
%% ============================================================

for k = 1:3

    filename = files{k};

    disp(['Reading: ' filename])

    %% Coordinates (read once)
    if k == 1
        lon   = double(ncread(filename,'longitude'));
        lat   = double(ncread(filename,'latitude'));
        depth = double(ncread(filename,'depth'));

        [Lon,Lat] = meshgrid(lon,lat);
    end

    %% Time
    time = double(ncread(filename,'time'));
    time = unixtime2mat(time);
    time = datetime(time,'ConvertFrom','datenum');

    %% NPP
    npp = double(ncread(filename,'nppv'));

    % lon x lat x depth x time

    %% Surface layer
    npp_surface = squeeze(npp(:,:,1,:));

    %% --------------------------------------------------------
    %% REMOVE CYCLONE DATES
    %% --------------------------------------------------------

    switch k

        case 1
            % Uesi + Gretel
            idx_keep = ~( ...
                (time >= datetime(2020,2,4)  & time <= datetime(2020,2,14)) | ...
                (time >= datetime(2020,3,14) & time <= datetime(2020,3,15)) );

        case 2
            % Kimi + Lucas + Niran
            idx_keep = ~( ...
                (time >= datetime(2021,1,16) & time <= datetime(2021,1,19)) | ...
                (time >= datetime(2021,1,24) & time <= datetime(2021,2,1))  | ...
                (time >= datetime(2021,2,27) & time <= datetime(2021,3,5)) );

        case 3
            % Ruby + Dovi
            idx_keep = ~( ...
                (time >= datetime(2021,12,9) & time <= datetime(2021,12,15)) | ...
                (time >= datetime(2022,2,7)  & time <= datetime(2022,2,12)) );
    end

    %% Clean data
    time_clean = time(idx_keep);
    npp_clean  = npp_surface(:,:,idx_keep);

    %% Concatenate all years
    baseline_all = cat(3,baseline_all,npp_clean);
    time_all     = [time_all; time_clean];

end

%% ============================================================
%% CHECK SIZE
%% ============================================================

disp('--- FINAL BASELINE SIZE ---')
size(baseline_all)

% lon x lat x total_days

%% ============================================================
%% MEAN AND STD
%% ============================================================

baseline_mean = mean(baseline_all,3,'omitnan');
baseline_std  = std(baseline_all,0,3,'omitnan');

baseline_mean = reshape(baseline_mean,length(lon),length(lat));
baseline_std  = reshape(baseline_std ,length(lon),length(lat));

%% ============================================================
%% TIMESERIES
%% ============================================================

ts_npp = squeeze(mean(mean(baseline_all,1,'omitnan'),2,'omitnan'));

%% ============================================================
%% MEAN MAP
%% ============================================================

figure;

m_proj('mercator',...
       'lon',[min(lon) max(lon)],...
       'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,baseline_mean')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12,'tickdir','out')

colorbar
colormap(jet)

title('Composite Baseline Mean Surface NPPV')

saveas(gcf,fullfile(output_dir,'BASELINE_NPPV_mean.png'))

%% ============================================================
%% STD MAP
%% ============================================================

figure;

m_proj('mercator',...
       'lon',[min(lon) max(lon)],...
       'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,baseline_std')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12,'tickdir','out')

colorbar
colormap(parula)

title('Composite Baseline STD Surface NPPV')

saveas(gcf,fullfile(output_dir,'BASELINE_NPPV_std.png'))

%% ============================================================
%% TIMESERIES
%% ============================================================

figure;

plot(time_all,ts_npp,'LineWidth',1.8)
grid on

xlabel('Date')
ylabel('NPPV')
title('Composite Baseline Daily Mean NPPV')

saveas(gcf,fullfile(output_dir,'BASELINE_NPPV_timeseries.png'))

%% ============================================================
%% PEAK EVENT MAP
%% ============================================================

[~,idx_peak] = max(ts_npp,[],'omitnan');

peak_date = time_all(idx_peak);
peak_map  = baseline_all(:,:,idx_peak);

figure;

m_proj('mercator',...
       'lon',[min(lon) max(lon)],...
       'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,peak_map')
shading interp
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12,'tickdir','out')

colorbar
colormap(jet)

title(['Composite Baseline Peak NPPV - ' datestr(peak_date)])

saveas(gcf,fullfile(output_dir,'BASELINE_NPPV_peakmap.png'))

%% ============================================================
%% LINEAR REGRESSION TREND
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
title('Composite Baseline Linear Trend')

legend('Observed','Linear fit','Location','best')

txt = sprintf('Slope = %.4f   R^2 = %.3f',p(1),R2);
text(min(x)+5,max(y)*0.95,txt,'FontSize',11)

saveas(gcf,fullfile(output_dir,'BASELINE_NPPV_regression.png'))

disp('--- REGRESSION RESULTS ---')
disp(['Slope = ' num2str(p(1))])
disp(['Intercept = ' num2str(p(2))])
disp(['R2 = ' num2str(R2)])