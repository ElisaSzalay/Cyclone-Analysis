%% Final Project: Analysis on Alfred and Fina cyclones compared to a baseline cyclone-free
% Variables analyzed: Chlorophyll-a (CHL) and nitrate (NO3)
% Variables dimensions: 
% CHL = lon x lat x time
% NO3 = lon x lat x depth x time
% NO3 averaged over 0–50 m
% Baseline years 
% Full area domain : 13-27°S, 150-165°N
% Overlap area : 18.5-19.5°S, 150-165°N

% Sections :
% 0. Settings
% 1. Time series in overlapping box
% 2. Anomaly time series
% 3. T-test for time series
% 4. Full domain maps
% 5. Anomaly maps:
% 1) 3-day mean anomaly around cyclone peak
% 2) 3-day significant anomaly, p < 0.05
% 6. Functions
%% 0. SETTINGS
clear; clc; close all;
set(groot,'defaultFigureColor','w'); %

% output folder
output_dir = fullfile(pwd,'plots');

% Overlap box
lat_box = [-19.5 -18.5];
lon_box = [155 156];

% no3 depth range
no3_depth_range = [0 50];

% Wind peaks dates
alfred_peak = datetime(2025,2,28);
fina_peak   = datetime(2011,12,21);

% FILES
alfred_ch  = 'Alfred_CH.nc';
alfred_no3 = 'Alfred_no3.nc';

fina_ch    = 'Fina_CH.nc';
fina_no3   = 'Fina_no3.nc';

baseline_ch_files  = {'baseline_CH_2019.nc','baseline_CH_2020.nc','baseline_CH_2021.nc'};
baseline_no3_files = {'baseline_no3_2019.nc','baseline_no3_2020.nc','baseline_no3_2021.nc'};

%% 1. TIME SERIES IN OVERLAP BOX

[ch_base_ts, t_base] = build_baseline(baseline_ch_files, ...
    @(f) extract_chl_box(f,lat_box,lon_box));

[no3_base_ts, ~] = build_baseline(baseline_no3_files, ...
    @(f) extract_no3_box_0_50m(f,lat_box,lon_box,no3_depth_range));

[ch_alfred_ts, t_alfred] = extract_chl_box(alfred_ch,lat_box,lon_box);
[no3_alfred_ts, t_alfred_no3] = extract_no3_box_0_50m(alfred_no3,lat_box,lon_box,no3_depth_range);

[ch_fina_ts, t_fina] = extract_chl_box(fina_ch,lat_box,lon_box);
[no3_fina_ts, t_fina_no3] = extract_no3_box_0_50m(fina_no3,lat_box,lon_box,no3_depth_range);

no3_alfred_ts = interp1(t_alfred_no3,no3_alfred_ts,t_alfred,'linear','extrap');
no3_fina_ts   = interp1(t_fina_no3,no3_fina_ts,t_fina,'linear','extrap');

plot_dual_timeseries(t_base, ch_base_ts, no3_base_ts, ...
    'Baseline period: November to March', ...
    'Baseline_CHL_NO3_0_50m', output_dir, NaT, false);

plot_dual_timeseries(t_alfred, ch_alfred_ts, no3_alfred_ts, ...
    'Cyclone Alfred: 01 Feb 2025 to 31 Mar 2025', ...
    'Alfred_CHL_NO3_0_50m', output_dir, alfred_peak, false);

plot_dual_timeseries(t_fina, ch_fina_ts, no3_fina_ts, ...
    'Cyclone Fina: 01 Dec 2011 to 31 Jan 2012', ...
    'Fina_CHL_NO3_0_50m', output_dir, fina_peak, false);

%% 2. ANOMALY TIME SERIES

[ch_base_for_alfred, no3_base_for_alfred] = baseline_for_event(t_alfred,t_base,ch_base_ts,no3_base_ts);
[ch_base_for_fina, no3_base_for_fina]     = baseline_for_event(t_fina,t_base,ch_base_ts,no3_base_ts);

anom_ch_alfred_ts  = ch_alfred_ts  - ch_base_for_alfred;
anom_no3_alfred_ts = no3_alfred_ts - no3_base_for_alfred;

anom_ch_fina_ts  = ch_fina_ts  - ch_base_for_fina;
anom_no3_fina_ts = no3_fina_ts - no3_base_for_fina;

plot_dual_timeseries(t_alfred, anom_ch_alfred_ts, anom_no3_alfred_ts, ...
    'Alfred anomaly relative to cyclone-free baseline', ...
    'Alfred_Anomaly_CHL_NO3_0_50m', output_dir, alfred_peak, true);

plot_dual_timeseries(t_fina, anom_ch_fina_ts, anom_no3_fina_ts, ...
    'Fina anomaly relative to cyclone-free baseline', ...
    'Fina_Anomaly_CHL_NO3_0_50m', output_dir, fina_peak, true);

%% 3. T-TESTS FOR TIME SERIES


[p_chl_alfred,t_chl_alfred] = manual_ttest2(ch_alfred_ts,ch_base_ts);
[p_no3_alfred,t_no3_alfred] = manual_ttest2(no3_alfred_ts,no3_base_ts);

[p_chl_fina,t_chl_fina] = manual_ttest2(ch_fina_ts,ch_base_ts);
[p_no3_fina,t_no3_fina] = manual_ttest2(no3_fina_ts,no3_base_ts);

[p_anom_chl_alfred,t_anom_chl_alfred] = manual_ttest2(anom_ch_alfred_ts,zeros(size(anom_ch_alfred_ts)));
[p_anom_no3_alfred,t_anom_no3_alfred] = manual_ttest2(anom_no3_alfred_ts,zeros(size(anom_no3_alfred_ts)));

[p_anom_chl_fina,t_anom_chl_fina] = manual_ttest2(anom_ch_fina_ts,zeros(size(anom_ch_fina_ts)));
[p_anom_no3_fina,t_anom_no3_fina] = manual_ttest2(anom_no3_fina_ts,zeros(size(anom_no3_fina_ts)));

disp(' ')
disp('============================================================')
disp('MANUAL T-TEST RESULTS')
disp('============================================================')

fprintf('\n--- ALFRED vs BASELINE ---\n')
fprintf('CHL: t = %.3f, p = %.5f -> %s\n', t_chl_alfred, p_chl_alfred, interpret_p(p_chl_alfred))
fprintf('NO3 0-50 m: t = %.3f, p = %.5f -> %s\n', t_no3_alfred, p_no3_alfred, interpret_p(p_no3_alfred))

fprintf('\n--- FINA vs BASELINE ---\n')
fprintf('CHL: t = %.3f, p = %.5f -> %s\n', t_chl_fina, p_chl_fina, interpret_p(p_chl_fina))
fprintf('NO3 0-50 m: t = %.3f, p = %.5f -> %s\n', t_no3_fina, p_no3_fina, interpret_p(p_no3_fina))

fprintf('\n--- ALFRED ANOMALY vs ZERO ---\n')
fprintf('CHL anomaly: t = %.3f, p = %.5f -> %s\n', t_anom_chl_alfred, p_anom_chl_alfred, interpret_p(p_anom_chl_alfred))
fprintf('NO3 0-50 m anomaly: t = %.3f, p = %.5f -> %s\n', t_anom_no3_alfred, p_anom_no3_alfred, interpret_p(p_anom_no3_alfred))

fprintf('\n--- FINA ANOMALY vs ZERO ---\n')
fprintf('CHL anomaly: t = %.3f, p = %.5f -> %s\n', t_anom_chl_fina, p_anom_chl_fina, interpret_p(p_anom_chl_fina))
fprintf('NO3 0-50 m anomaly: t = %.3f, p = %.5f -> %s\n', t_anom_no3_fina, p_anom_no3_fina, interpret_p(p_anom_no3_fina))

ttest_results = table( ...
    ["Alfred vs Baseline"; "Alfred vs Baseline"; ...
     "Fina vs Baseline"; "Fina vs Baseline"; ...
     "Alfred anomaly vs zero"; "Alfred anomaly vs zero"; ...
     "Fina anomaly vs zero"; "Fina anomaly vs zero"], ...
    ["CHL"; "NO3_0_50m"; "CHL"; "NO3_0_50m"; ...
     "CHL anomaly"; "NO3_0_50m anomaly"; "CHL anomaly"; "NO3_0_50m anomaly"], ...
    [t_chl_alfred; t_no3_alfred; t_chl_fina; t_no3_fina; ...
     t_anom_chl_alfred; t_anom_no3_alfred; t_anom_chl_fina; t_anom_no3_fina], ...
    [p_chl_alfred; p_no3_alfred; p_chl_fina; p_no3_fina; ...
     p_anom_chl_alfred; p_anom_no3_alfred; p_anom_chl_fina; p_anom_no3_fina], ...
    [p_chl_alfred<0.05; p_no3_alfred<0.05; ...
     p_chl_fina<0.05; p_no3_fina<0.05; ...
     p_anom_chl_alfred<0.05; p_anom_no3_alfred<0.05; ...
     p_anom_chl_fina<0.05; p_anom_no3_fina<0.05], ...
    'VariableNames',{'Test','Variable','t_value','p_value','Significant_005'} );

writetable(ttest_results,fullfile(output_dir,'ttest_results_NO3_0_50m.csv'));


%% 4. FULL-DOMAIN MAPS
lat_ch = ncread(alfred_ch,'latitude');
lon_ch = ncread(alfred_ch,'longitude');

lat_no3 = ncread(alfred_no3,'latitude');
lon_no3 = ncread(alfred_no3,'longitude');

% Baseline maps
[ch_base_mean,ch_base_std] = baseline_maps_chl(baseline_ch_files);
[no3_base_mean,no3_base_std] = baseline_maps_no3_0_50m(baseline_no3_files,no3_depth_range);

% Alfred mean/std maps
[ch_mean_a,ch_std_a] = spatial_stats_chl(alfred_ch);
[no3_mean_a,no3_std_a] = spatial_stats_no3_0_50m(alfred_no3,no3_depth_range);

% Fina mean/std maps
[ch_mean_f,ch_std_f] = spatial_stats_chl(fina_ch);
[no3_mean_f,no3_std_f] = spatial_stats_no3_0_50m(fina_no3,no3_depth_range);


%% 5. 3-DAY MEAN ANOMALY + GRID-CELL SIGNIFICANCE TEST

% Alfred 3-day cyclone mean
ch_peak3_a = peak3_map_chl(alfred_ch,alfred_peak);
no3_peak3_a = peak3_map_no3_0_50m(alfred_no3,alfred_peak,no3_depth_range);

% Alfred 3-day cyclone-free baseline stack
ch_base_stack_a = baseline_3day_stack_chl(baseline_ch_files,alfred_peak);
no3_base_stack_a = baseline_3day_stack_no3_0_50m(baseline_no3_files,alfred_peak,no3_depth_range);

% Alfred 3-day baseline mean
ch_base_peak3_a = mean(ch_base_stack_a,3,'omitnan');
no3_base_peak3_a = mean(no3_base_stack_a,3,'omitnan');

% Alfred 3-day anomalies
ch_peak3_anom_a  = ch_peak3_a  - ch_base_peak3_a;
no3_peak3_anom_a = no3_peak3_a - no3_base_peak3_a;

% Alfred 3-day significant anomalies
[ch_peak3_sig_a, ch_pmap_a] = gridcell_ttest_mask(ch_peak3_anom_a,ch_peak3_a,ch_base_stack_a,0.05);
[no3_peak3_sig_a, no3_pmap_a] = gridcell_ttest_mask(no3_peak3_anom_a,no3_peak3_a,no3_base_stack_a,0.05);

% Fina 3-day cyclone mean
ch_peak3_f = peak3_map_chl(fina_ch,fina_peak);
no3_peak3_f = peak3_map_no3_0_50m(fina_no3,fina_peak,no3_depth_range);

% Fina 3-day cyclone-free baseline stack
ch_base_stack_f = baseline_3day_stack_chl(baseline_ch_files,fina_peak);
no3_base_stack_f = baseline_3day_stack_no3_0_50m(baseline_no3_files,fina_peak,no3_depth_range);

% Fina 3-day baseline mean
ch_base_peak3_f = mean(ch_base_stack_f,3,'omitnan');
no3_base_peak3_f = mean(no3_base_stack_f,3,'omitnan');

% Fina 3-day anomalies
ch_peak3_anom_f  = ch_peak3_f  - ch_base_peak3_f;
no3_peak3_anom_f = no3_peak3_f - no3_base_peak3_f;

% Fina 3-day significant anomalies
[ch_peak3_sig_f, ch_pmap_f] = gridcell_ttest_mask(ch_peak3_anom_f,ch_peak3_f,ch_base_stack_f,0.05);
[no3_peak3_sig_f, no3_pmap_f] = gridcell_ttest_mask(no3_peak3_anom_f,no3_peak3_f,no3_base_stack_f,0.05);

save(fullfile(output_dir,'gridcell_pmaps.mat'), ...
    'ch_pmap_a','no3_pmap_a','ch_pmap_f','no3_pmap_f');

%% Color limits
clims.ch_mean = safe_clim([0 max([ch_mean_a(:); ch_mean_f(:)],[],'omitnan')]);
clims.no3_mean = safe_clim([0 max([no3_mean_a(:); no3_mean_f(:)],[],'omitnan')]);

clims.ch_std = safe_clim([0 max([ch_std_a(:); ch_std_f(:)],[],'omitnan')]);
clims.no3_std = safe_clim([0 max([no3_std_a(:); no3_std_f(:)],[],'omitnan')]);

ch_anom_abs = percentile_no_toolbox(abs([ch_peak3_anom_a(:); ch_peak3_anom_f(:); ch_peak3_sig_a(:); ch_peak3_sig_f(:)]),98);
no3_anom_abs = percentile_no_toolbox(abs([no3_peak3_anom_a(:); no3_peak3_anom_f(:); no3_peak3_sig_a(:); no3_peak3_sig_f(:)]),98);

clims.ch_anom = safe_clim([-ch_anom_abs ch_anom_abs]);
clims.no3_anom = safe_clim([-no3_anom_abs no3_anom_abs]);

%% Mean and STD figures
plot_mean_std_comparison(ch_mean_a,ch_std_a,ch_mean_f,ch_std_f, ...
    lat_ch,lon_ch,lat_box,lon_box, ...
    'CHL Mean and Standard Deviation', ...
    'CHL_Mean_STD_Comparison', ...
    'mg m^{-3}', output_dir, clims.ch_mean, clims.ch_std);

plot_mean_std_comparison(no3_mean_a,no3_std_a,no3_mean_f,no3_std_f, ...
    lat_no3,lon_no3,lat_box,lon_box, ...
    'NO_3 0-50 m Mean and Standard Deviation', ...
    'NO3_0_50m_Mean_STD_Comparison', ...
    'mmol m^{-3}', output_dir, clims.no3_mean, clims.no3_std);

%% 3-day anomaly figures
plot_anomaly_comparison(ch_peak3_anom_a,ch_peak3_sig_a,ch_peak3_anom_f,ch_peak3_sig_f, ...
    lat_ch,lon_ch,lat_box,lon_box, ...
    'CHL 3-day Anomaly Comparison', ...
    'CHL_3day_Anomaly_Comparison', ...
    'mg m^{-3}', output_dir, clims.ch_anom);

plot_anomaly_comparison(no3_peak3_anom_a,no3_peak3_sig_a,no3_peak3_anom_f,no3_peak3_sig_f, ...
    lat_no3,lon_no3,lat_box,lon_box, ...
    'NO_3 0-50 m 3-day Anomaly Comparison', ...
    'NO3_0_50m_3day_Anomaly_Comparison', ...
    'mmol m^{-3}', output_dir, clims.no3_anom);

disp('Analysis completed successfully.');
disp("All figures and t-test table exported to: " + output_dir);


%% 6. Functions
function [series,time] = extract_chl_box(file,lat_box,lon_box)

lat = ncread(file,'latitude');
lon = ncread(file,'longitude');
chl = ncread(file,'CHL');
time = datetime(double(ncread(file,'time')),'ConvertFrom','posixtime');

lon_i = lon >= lon_box(1) & lon <= lon_box(2);
lat_i = lat >= lat_box(1) & lat <= lat_box(2);

nt = size(chl,3);
series = NaN(nt,1);

for t = 1:nt
    frame = chl(:,:,t);
    sub = frame(lon_i,lat_i);
    series(t) = mean(sub,'all','omitnan');
end

series = series(:);
time = time(:);

end

function [series,time] = extract_no3_box_0_50m(file,lat_box,lon_box,depth_range)

lat = ncread(file,'latitude');
lon = ncread(file,'longitude');
depth = ncread(file,'depth');
no3 = ncread(file,'no3');
time = datetime(double(ncread(file,'time')),'ConvertFrom','posixtime');

lon_i = lon >= lon_box(1) & lon <= lon_box(2);
lat_i = lat >= lat_box(1) & lat <= lat_box(2);
dep_i = depth >= depth_range(1) & depth <= depth_range(2);

nt = size(no3,4);
series = NaN(nt,1);

for t = 1:nt
    frame = no3(:,:,:,t);
    sub = frame(lon_i,lat_i,dep_i);
    series(t) = mean(sub,'all','omitnan');
end

series = series(:);
time = time(:);

end

function [baseline,t_base] = build_baseline(files,extractor)

exclude_periods = get_exclude_periods();

t_base = (datetime(2001,11,1):days(1):datetime(2002,3,31))';
nDays = length(t_base);

all_years_data = NaN(nDays,length(files));

for i = 1:length(files)

    [s,t] = extractor(files{i});
    s = s(:);
    t = t(:);

    keep = month(t) >= 11 | month(t) <= 3;
    s = s(keep);
    t = t(keep);

    remove_mask = false(size(t));

    for j = 1:size(exclude_periods,1)
        remove_mask = remove_mask | ...
            (t >= exclude_periods(j,1) & t <= exclude_periods(j,2));
    end

    t(remove_mask) = [];
    s(remove_mask) = [];

    for k = 1:length(t)

        if month(t(k)) >= 11
            ref_date = datetime(2001,month(t(k)),day(t(k)));
        else
            ref_date = datetime(2002,month(t(k)),day(t(k)));
        end

        season_day = days(ref_date - datetime(2001,11,1)) + 1;

        if season_day >= 1 && season_day <= nDays
            all_years_data(season_day,i) = s(k);
        end

    end
end

baseline = mean(all_years_data,2,'omitnan');

end

function [base_ch,base_no3] = baseline_for_event(t_event,t_base,ch_base,no3_base)

base_ch = NaN(length(t_event),1);
base_no3 = NaN(length(t_event),1);

for i = 1:length(t_event)

    m = month(t_event(i));
    d = day(t_event(i));

    if m >= 11
        ref_date = datetime(2001,m,d);
    else
        ref_date = datetime(2002,m,d);
    end

    idx = find(t_base == ref_date,1);

    if ~isempty(idx)
        base_ch(i) = ch_base(idx);
        base_no3(i) = no3_base(idx);
    end

end

end

function plot_dual_timeseries(time,chl,no3,title_str,filename,output_dir,peak_date,show_zero_line)

figure;

yyaxis left
h1 = plot(time,chl,'LineWidth',1.5);
ylabel('Chlorophyll-a (mg m^{-3})')
hold on

if show_zero_line
    yline(0,'--r','LineWidth',1.2);
end

yyaxis right
h2 = plot(time,no3,'LineWidth',1.5);
ylabel('Nitrate NO_3 0-50 m (mmol m^{-3})')
hold on

if show_zero_line
    yline(0,'--r','LineWidth',1.2);
end

if ~isnat(peak_date)
    xline(peak_date,'--','Color',[0.45 0 0],'LineWidth',1.5, ...
        'Label','Wind peak','LabelOrientation','horizontal', ...
        'LabelVerticalAlignment','bottom');
end

xlabel('Time')
title("CHL and NO_3 0-50 m during " + title_str)
legend([h1 h2],{'CHL','NO_3 0-50 m'},'Location','best')
grid on
box on
set(gca,'FontSize',10)
xtickformat('dd-MMM')

set(gcf,'Toolbar','none')
set(gcf,'MenuBar','none')

exportgraphics(gcf,fullfile(output_dir,filename + ".png"),'Resolution',300);
exportgraphics(gcf,fullfile(output_dir,filename + ".pdf"),'ContentType','vector');

end

function [mean_map,std_map] = spatial_stats_chl(file)

data = ncread(file,'CHL');
mean_map = mean(data,3,'omitnan');
std_map  = std(data,0,3,'omitnan');

end

function [mean_map,std_map] = spatial_stats_no3_0_50m(file,depth_range)

data = ncread(file,'no3');
depth = ncread(file,'depth');

dep_i = depth >= depth_range(1) & depth <= depth_range(2);

data_0_50 = data(:,:,dep_i,:);
data_0_50 = squeeze(mean(data_0_50,3,'omitnan'));

mean_map = mean(data_0_50,3,'omitnan');
std_map  = std(data_0_50,0,3,'omitnan');

end

function [baseline_mean,baseline_std] = baseline_maps_chl(files)

exclude_periods = get_exclude_periods();
all_maps = [];

for i = 1:length(files)

    data = ncread(files{i},'CHL');
    time = datetime(double(ncread(files{i},'time')),'ConvertFrom','posixtime');

    keep = month(time) >= 11 | month(time) <= 3;

    for j = 1:size(exclude_periods,1)
        keep = keep & ~(time >= exclude_periods(j,1) & time <= exclude_periods(j,2));
    end

    data = data(:,:,keep);
    all_maps = cat(3,all_maps,data);

end

baseline_mean = mean(all_maps,3,'omitnan');
baseline_std  = std(all_maps,0,3,'omitnan');

end

function [baseline_mean,baseline_std] = baseline_maps_no3_0_50m(files,depth_range)

exclude_periods = get_exclude_periods();
all_maps = [];

for i = 1:length(files)

    data = ncread(files{i},'no3');
    depth = ncread(files{i},'depth');
    time = datetime(double(ncread(files{i},'time')),'ConvertFrom','posixtime');

    keep = month(time) >= 11 | month(time) <= 3;

    for j = 1:size(exclude_periods,1)
        keep = keep & ~(time >= exclude_periods(j,1) & time <= exclude_periods(j,2));
    end

    dep_i = depth >= depth_range(1) & depth <= depth_range(2);

    data_0_50 = data(:,:,dep_i,keep);
    data_0_50 = squeeze(mean(data_0_50,3,'omitnan'));

    all_maps = cat(3,all_maps,data_0_50);

end

baseline_mean = mean(all_maps,3,'omitnan');
baseline_std  = std(all_maps,0,3,'omitnan');

end

function peak3 = peak3_map_chl(file,peak_date)

data = ncread(file,'CHL');
time = datetime(double(ncread(file,'time')),'ConvertFrom','posixtime');

idx = time >= peak_date-days(1) & time <= peak_date+days(1);
peak3 = mean(data(:,:,idx),3,'omitnan');

end

function peak3 = peak3_map_no3_0_50m(file,peak_date,depth_range)

data = ncread(file,'no3');
depth = ncread(file,'depth');
time = datetime(double(ncread(file,'time')),'ConvertFrom','posixtime');

dep_i = depth >= depth_range(1) & depth <= depth_range(2);
idx = time >= peak_date-days(1) & time <= peak_date+days(1);

tmp = data(:,:,dep_i,idx);
tmp = squeeze(mean(tmp,3,'omitnan'));

peak3 = mean(tmp,3,'omitnan');

end

function baseline_stack = baseline_3day_stack_chl(files,event_date)

exclude_periods = get_exclude_periods();
baseline_stack = [];

target_dates = event_date + days(-1:1);

for i = 1:length(files)

    data = ncread(files{i},'CHL');
    time = datetime(double(ncread(files{i},'time')),'ConvertFrom','posixtime');

    keep = false(size(time));

    for d = 1:length(target_dates)
        keep = keep | (month(time)==month(target_dates(d)) & day(time)==day(target_dates(d)));
    end

    for j = 1:size(exclude_periods,1)
        keep = keep & ~(time >= exclude_periods(j,1) & time <= exclude_periods(j,2));
    end

    if any(keep)
        one_season_mean = mean(data(:,:,keep),3,'omitnan');
        baseline_stack = cat(3,baseline_stack,one_season_mean);
    end

end

end

function baseline_stack = baseline_3day_stack_no3_0_50m(files,event_date,depth_range)

exclude_periods = get_exclude_periods();
baseline_stack = [];

target_dates = event_date + days(-1:1);

for i = 1:length(files)

    data = ncread(files{i},'no3');
    depth = ncread(files{i},'depth');
    time = datetime(double(ncread(files{i},'time')),'ConvertFrom','posixtime');

    dep_i = depth >= depth_range(1) & depth <= depth_range(2);

    keep = false(size(time));

    for d = 1:length(target_dates)
        keep = keep | (month(time)==month(target_dates(d)) & day(time)==day(target_dates(d)));
    end

    for j = 1:size(exclude_periods,1)
        keep = keep & ~(time >= exclude_periods(j,1) & time <= exclude_periods(j,2));
    end

    if any(keep)
        tmp = data(:,:,dep_i,keep);
        tmp = squeeze(mean(tmp,3,'omitnan'));
        one_season_mean = mean(tmp,3,'omitnan');
        baseline_stack = cat(3,baseline_stack,one_season_mean);
    end

end

end

function [sig_anom,pmap] = gridcell_ttest_mask(anom_map,cyclone_map,baseline_stack,alpha)

baseline_mean = mean(baseline_stack,3,'omitnan');
baseline_std  = std(baseline_stack,0,3,'omitnan');

n = sum(~isnan(baseline_stack),3);
se = baseline_std ./ sqrt(n);

tmap = (cyclone_map - baseline_mean) ./ se;

pmap = erfc(abs(tmap)./sqrt(2));

sig_mask = pmap < alpha;
sig_mask(isnan(pmap)) = false;
sig_mask(n < 2) = false;
sig_mask(se == 0) = false;

sig_anom = anom_map;
sig_anom(~sig_mask) = NaN;

end

function plot_mean_std_comparison(mean_a,std_a,mean_f,std_f,lat,lon,lat_box,lon_box,title_str,filename,var_label,output_dir,mean_lim,std_lim)

figure('Position',[100 100 1100 800]);
set(gcf,'Toolbar','none')
set(gcf,'MenuBar','none')

land_color = [0.78 0.62 0.42];
cmap = turbo(256);

tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile
plot_single_map(mean_a,lat,lon,lat_box,lon_box,cmap,land_color,mean_lim)
title('Alfred Mean')
cb = colorbar; cb.Label.String = var_label;

nexttile
plot_single_map(std_a,lat,lon,lat_box,lon_box,cmap,land_color,std_lim)
title('Alfred Standard Deviation')
cb = colorbar; cb.Label.String = var_label;

nexttile
plot_single_map(mean_f,lat,lon,lat_box,lon_box,cmap,land_color,mean_lim)
title('Fina Mean')
cb = colorbar; cb.Label.String = var_label;

nexttile
plot_single_map(std_f,lat,lon,lat_box,lon_box,cmap,land_color,std_lim)
title('Fina Standard Deviation')
cb = colorbar; cb.Label.String = var_label;

sgtitle(title_str)

exportgraphics(gcf,fullfile(output_dir,filename + ".png"),'Resolution',300);
exportgraphics(gcf,fullfile(output_dir,filename + ".pdf"),'ContentType','vector');

end

function plot_anomaly_comparison(anom_a,sig_a,anom_f,sig_f,lat,lon,lat_box,lon_box,title_str,filename,var_label,output_dir,anom_lim)

figure('Position',[100 100 1100 800]);
set(gcf,'Toolbar','none')
set(gcf,'MenuBar','none')

land_color = [0.78 0.62 0.42];
nonsig_color = [0.96 0.96 0.96];
cmap = redblue(256);

tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile
plot_single_map(anom_a,lat,lon,lat_box,lon_box,cmap,land_color,anom_lim)
title('Alfred 3-day mean anomaly')
cb = colorbar; cb.Label.String = var_label;

nexttile
plot_single_map(sig_a,lat,lon,lat_box,lon_box,cmap,nonsig_color,anom_lim)
title('Alfred 3-day significant anomaly')
cb = colorbar; cb.Label.String = var_label;

nexttile
plot_single_map(anom_f,lat,lon,lat_box,lon_box,cmap,land_color,anom_lim)
title('Fina 3-day mean anomaly')
cb = colorbar; cb.Label.String = var_label;

nexttile
plot_single_map(sig_f,lat,lon,lat_box,lon_box,cmap,nonsig_color,anom_lim)
title('Fina 3-day significant anomaly')
cb = colorbar; cb.Label.String = var_label;

sgtitle(title_str)

exportgraphics(gcf,fullfile(output_dir,filename + ".png"),'Resolution',300);
exportgraphics(gcf,fullfile(output_dir,filename + ".pdf"),'ContentType','vector');

end

function plot_single_map(map_data,lat,lon,lat_box,lon_box,map_cmap,background_color,clim_values)

h = imagesc(lon,lat,map_data');
set(gca,'YDir','normal')
set(gca,'Color',background_color)
set(h,'AlphaData',~isnan(map_data'))

colormap(gca,map_cmap)
caxis(clim_values)

xlabel('Longitude (°E)')
ylabel('Latitude (°N)')
grid on
box on
hold on

x_box = [lon_box(1) lon_box(2) lon_box(2) lon_box(1) lon_box(1)];
y_box = [lat_box(1) lat_box(1) lat_box(2) lat_box(2) lat_box(1)];

plot(x_box,y_box,'r-','LineWidth',2)

end

function [p,t] = manual_ttest2(x,y)

x = x(~isnan(x));
y = y(~isnan(y));

nx = length(x);
ny = length(y);

if nx < 2 || ny < 2
    p = NaN;
    t = NaN;
    return
end

mx = mean(x);
my = mean(y);

vx = var(x);
vy = var(y);

se = sqrt(vx/nx + vy/ny);

if se == 0 || isnan(se)
    p = NaN;
    t = NaN;
    return
end

t = (mx - my) / se;
p = erfc(abs(t)/sqrt(2));

end

function txt = interpret_p(p)

if isnan(p)
    txt = 'not available';
elseif p < 0.01
    txt = 'highly significant';
elseif p < 0.05
    txt = 'significant';
else
    txt = 'not significant';
end

end

function exclude_periods = get_exclude_periods()

exclude_periods = [
    datetime(2022,2,7)   datetime(2022,2,12);
    datetime(2021,12,9)  datetime(2021,12,15);
    datetime(2021,2,27)  datetime(2021,3,5);
    datetime(2021,1,24)  datetime(2021,2,1);
    datetime(2021,1,16)  datetime(2021,1,19);
    datetime(2020,4,1)   datetime(2020,4,10);
    datetime(2020,3,14)  datetime(2020,3,15);
    datetime(2020,2,4)   datetime(2020,2,14);
];

end

function clim_out = safe_clim(clim_in)

clim_out = clim_in;

if any(isnan(clim_out)) || any(isinf(clim_out))
    clim_out = [0 1];
end

if clim_out(1) == clim_out(2)
    clim_out = clim_out + [-1 1]*eps;
end

end

function p = percentile_no_toolbox(x,perc)

x = x(~isnan(x));
x = sort(x(:));

if isempty(x)
    p = NaN;
    return
end

idx = round((perc/100) * length(x));
idx = max(1,min(idx,length(x)));

p = x(idx);

end

function cmap = redblue(n)

if nargin < 1
    n = 256;
end

bottom = [0 0 0.5];
middle = [1 1 1];
top = [0.5 0 0];

n1 = floor(n/2);
n2 = n - n1;

cmap = [
    linspace(bottom(1),middle(1),n1)' linspace(bottom(2),middle(2),n1)' linspace(bottom(3),middle(3),n1)';
    linspace(middle(1),top(1),n2)' linspace(middle(2),top(2),n2)' linspace(middle(3),top(3),n2)'
];

end