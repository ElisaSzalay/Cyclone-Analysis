%% =========================================================================
%  UNIFIED ANALYSIS — NPPV & SPM
%  Cyclone Alfred vs Cyclone Fina vs Baseline
%  -------------------------------------------------------------------------
%  Variables:
%    - NPPV : Net Primary Production  [lon x lat x depth x time]  (surface)
%    - SPM  : Suspended Particulate Matter  [lon x lat x time]
%
%
%  Structure:
%    0. Settings
%    1. Load baseline data
%    2. Load Alfred & Fina data
%    3. Compute statistics (mean, std, anomaly, z-score)
%    4. T-test with 3-day means around cyclone peak
%    5. Figures 
%    6. Save results
%    7. Local functions
% =========================================================================

clear; clc; close all;

%% =========================================================================
%  0. SETTINGS
% =========================================================================

% ---- Cyclone peak dates (date of maximum wind intensity) ----
alfred_peak = datetime(2025, 2, 28, 'TimeZone', 'UTC');
fina_peak   = datetime(2011, 12, 21, 'TimeZone', 'UTC');

% ---- Cyclone time windows ----
alfred_start = datetime(2025, 2,  1, 'TimeZone', 'UTC');
alfred_end   = datetime(2025, 3, 31, 'TimeZone', 'UTC');
fina_start   = datetime(2011, 12,  1, 'TimeZone', 'UTC');
fina_end     = datetime(2012,  1, 31, 'TimeZone', 'UTC');

% ---- Overlap analysis box ----
lat_box = [-19.5, -18.5];
lon_box = [155.0, 156.0];

% ---- NetCDF file names ----
alfred_nppv_file = 'nppv_alfred_bo.nc';
alfred_spm_file  = 'spm_alfred_bo.nc';
fina_nppv_file   = 'nppv_fina_bo.nc';
fina_spm_file    = 'spm_fina_bo.nc';

baseline_nppv_files = {'nppv_baseline_bo1.nc', 'nppv_baseline_bo2.nc', 'nppv_baseline_bo3.nc'};
baseline_spm_files  = {'spm_baseline_bo1.nc',  'spm_baseline_bo2.nc',  'spm_baseline_bo3.nc'};

% ---- Variable names inside the NetCDF files ----
nppv_var = 'nppv';
spm_var  = 'SPM';    % uppercase — confirmed from nc files
lon_var  = 'longitude';
lat_var  = 'latitude';
time_var = 'time';

% ---- Output folder ----
output_dir = fullfile(pwd, 'plots_NPPV_SPM');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% ---- Figure style ----
set(groot, 'defaultFigureColor', 'w');
LAND_COLOR   = [0.78 0.62 0.42];  % beige for land / background
NONSIG_COLOR = [0.92 0.92 0.92];  % grey for non-significant pixels

fprintf('=== UNIFIED NPPV & SPM ANALYSIS ===\n\n');

%% =========================================================================
%  1. LOAD BASELINE DATA
%  FIX (E): Cyclone exclusions are now applied inside load_baseline_data
%  via get_exclude_periods(), which contains all 8 forbidden periods.
% =========================================================================

fprintf('--- Loading baseline NPPV data ---\n');
[all_nppv_base, all_time_nppv_base, lon_nppv, lat_nppv] = ...
    load_baseline_data(baseline_nppv_files, nppv_var, lon_var, lat_var, time_var, true);

fprintf('--- Loading baseline SPM data ---\n');
[all_spm_base, all_time_spm_base, lon_spm, lat_spm] = ...
    load_baseline_data(baseline_spm_files, spm_var, lon_var, lat_var, time_var, false);

fprintf('Baseline loaded: NPPV = %d timesteps | SPM = %d timesteps\n\n', ...
    size(all_nppv_base, 3), size(all_spm_base, 3));

% Baseline spatial statistics
base_nppv_mean = mean(all_nppv_base, 3, 'omitnan');
base_nppv_std  = std(all_nppv_base,  0, 3, 'omitnan');
base_spm_mean  = mean(all_spm_base,  3, 'omitnan');
base_spm_std   = std(all_spm_base,   0, 3, 'omitnan');

% Baseline box time series (spatial mean within overlap box)
ts_base_nppv = extract_box_ts(all_nppv_base, lon_nppv, lat_nppv, lon_box, lat_box);
ts_base_spm  = extract_box_ts(all_spm_base,  lon_spm,  lat_spm,  lon_box, lat_box);

%% =========================================================================
%  2. LOAD ALFRED & FINA DATA
% =========================================================================

fprintf('--- Loading Cyclone Alfred data ---\n');
[nppv_alfred, t_nppv_alfred] = load_cyclone_data(alfred_nppv_file, nppv_var, ...
    time_var, alfred_start, alfred_end, true);
[spm_alfred,  t_spm_alfred]  = load_cyclone_data(alfred_spm_file,  spm_var,  ...
    time_var, alfred_start, alfred_end, false);

fprintf('Alfred loaded: NPPV = %d timesteps | SPM = %d timesteps\n\n', ...
    size(nppv_alfred, 3), size(spm_alfred, 3));

fprintf('--- Loading Cyclone Fina data ---\n');
[nppv_fina, t_nppv_fina] = load_cyclone_data(fina_nppv_file, nppv_var, ...
    time_var, fina_start, fina_end, true);
[spm_fina,  t_spm_fina]  = load_cyclone_data(fina_spm_file,  spm_var,  ...
    time_var, fina_start, fina_end, false);

fprintf('Fina loaded: NPPV = %d timesteps | SPM = %d timesteps\n\n', ...
    size(nppv_fina, 3), size(spm_fina, 3));

%% =========================================================================
%  3. COMPUTE STATISTICS
% =========================================================================

fprintf('--- Computing statistics ---\n');

% Spatial mean and std maps for each cyclone
alfred_nppv_mean = mean(nppv_alfred, 3, 'omitnan');
alfred_nppv_std  = std(nppv_alfred,  0, 3, 'omitnan');
alfred_spm_mean  = mean(spm_alfred,  3, 'omitnan');
alfred_spm_std   = std(spm_alfred,   0, 3, 'omitnan');

fina_nppv_mean   = mean(nppv_fina,   3, 'omitnan');
fina_nppv_std    = std(nppv_fina,    0, 3, 'omitnan');
fina_spm_mean    = mean(spm_fina,    3, 'omitnan');
fina_spm_std     = std(spm_fina,     0, 3, 'omitnan');

% Anomaly maps: cyclone_mean - baseline_mean
anom_nppv_alfred = alfred_nppv_mean - base_nppv_mean;
anom_nppv_fina   = fina_nppv_mean   - base_nppv_mean;
anom_spm_alfred  = alfred_spm_mean  - base_spm_mean;
anom_spm_fina    = fina_spm_mean    - base_spm_mean;

% Z-score maps: anomaly / baseline_std
zscore_nppv_alfred = anom_nppv_alfred ./ base_nppv_std;
zscore_nppv_fina   = anom_nppv_fina   ./ base_nppv_std;
zscore_spm_alfred  = anom_spm_alfred  ./ base_spm_std;
zscore_spm_fina    = anom_spm_fina    ./ base_spm_std;

% Box time series for each cyclone (spatial mean in overlap box)
ts_nppv_alfred = extract_box_ts(nppv_alfred, lon_nppv, lat_nppv, lon_box, lat_box);
ts_nppv_fina   = extract_box_ts(nppv_fina,   lon_nppv, lat_nppv, lon_box, lat_box);
ts_spm_alfred  = extract_box_ts(spm_alfred,  lon_spm,  lat_spm,  lon_box, lat_box);
ts_spm_fina    = extract_box_ts(spm_fina,    lon_spm,  lat_spm,  lon_box, lat_box);

% Anomaly time series: match baseline value to each cyclone date by season day
ts_base_nppv_for_alfred = interp_baseline_to_event(t_nppv_alfred, all_time_nppv_base, ts_base_nppv);
ts_base_nppv_for_fina   = interp_baseline_to_event(t_nppv_fina,   all_time_nppv_base, ts_base_nppv);
ts_base_spm_for_alfred  = interp_baseline_to_event(t_spm_alfred,  all_time_spm_base,  ts_base_spm);
ts_base_spm_for_fina    = interp_baseline_to_event(t_spm_fina,    all_time_spm_base,  ts_base_spm);

anom_ts_nppv_alfred = ts_nppv_alfred - ts_base_nppv_for_alfred;
anom_ts_nppv_fina   = ts_nppv_fina   - ts_base_nppv_for_fina;
anom_ts_spm_alfred  = ts_spm_alfred  - ts_base_spm_for_alfred;
anom_ts_spm_fina    = ts_spm_fina    - ts_base_spm_for_fina;

fprintf('  Alfred NPPV anomaly (domain mean): %.4f\n', mean(anom_nppv_alfred(:), 'omitnan'));
fprintf('  Alfred SPM  anomaly (domain mean): %.4f\n', mean(anom_spm_alfred(:),  'omitnan'));
fprintf('  Fina   NPPV anomaly (domain mean): %.4f\n', mean(anom_nppv_fina(:),   'omitnan'));
fprintf('  Fina   SPM  anomaly (domain mean): %.4f\n', mean(anom_spm_fina(:),    'omitnan'));

%% =========================================================================
%  4. T-TEST WITH 3-DAY MEANS AROUND CYCLONE PEAK


fprintf('\n--- T-TEST: 3-day peak means vs baseline (grid cell level) ---\n');

% Build 3-day cyclone peak maps
peak3_nppv_alfred = build_3day_peak_map(nppv_alfred, t_nppv_alfred, alfred_peak);
peak3_spm_alfred  = build_3day_peak_map(spm_alfred,  t_spm_alfred,  alfred_peak);
peak3_nppv_fina   = build_3day_peak_map(nppv_fina,   t_nppv_fina,   fina_peak);
peak3_spm_fina    = build_3day_peak_map(spm_fina,    t_spm_fina,    fina_peak);

% Build baseline 3-day stacks (same calendar days, all years)
base_stack_nppv_alfred = build_baseline_3day_stack(all_nppv_base, all_time_nppv_base, alfred_peak);
base_stack_spm_alfred  = build_baseline_3day_stack(all_spm_base,  all_time_spm_base,  alfred_peak);
base_stack_nppv_fina   = build_baseline_3day_stack(all_nppv_base, all_time_nppv_base, fina_peak);
base_stack_spm_fina    = build_baseline_3day_stack(all_spm_base,  all_time_spm_base,  fina_peak);

% 3-day anomaly maps
peak3_anom_nppv_alfred = peak3_nppv_alfred - mean(base_stack_nppv_alfred, 3, 'omitnan');
peak3_anom_spm_alfred  = peak3_spm_alfred  - mean(base_stack_spm_alfred,  3, 'omitnan');
peak3_anom_nppv_fina   = peak3_nppv_fina   - mean(base_stack_nppv_fina,   3, 'omitnan');
peak3_anom_spm_fina    = peak3_spm_fina    - mean(base_stack_spm_fina,    3, 'omitnan');

% Grid-cell t-test: significant anomaly mask
[sig_nppv_alfred, pmap_nppv_alfred] = gridcell_ttest_mask(peak3_nppv_alfred, base_stack_nppv_alfred, 0.05);
[sig_spm_alfred,  pmap_spm_alfred]  = gridcell_ttest_mask(peak3_spm_alfred,  base_stack_spm_alfred,  0.05);
[sig_nppv_fina,   pmap_nppv_fina]   = gridcell_ttest_mask(peak3_nppv_fina,   base_stack_nppv_fina,   0.05);
[sig_spm_fina,    pmap_spm_fina]    = gridcell_ttest_mask(peak3_spm_fina,    base_stack_spm_fina,    0.05);

% Apply sig mask to anomaly (NaN where not significant)
sig_anom_nppv_alfred = peak3_anom_nppv_alfred; sig_anom_nppv_alfred(isnan(sig_nppv_alfred)) = NaN;
sig_anom_spm_alfred  = peak3_anom_spm_alfred;  sig_anom_spm_alfred(isnan(sig_spm_alfred))   = NaN;
sig_anom_nppv_fina   = peak3_anom_nppv_fina;   sig_anom_nppv_fina(isnan(sig_nppv_fina))     = NaN;
sig_anom_spm_fina    = peak3_anom_spm_fina;    sig_anom_spm_fina(isnan(sig_spm_fina))       = NaN;

% ---- Domain-box t-test (time series) ----
fprintf('\n--- T-TEST: Domain box time series (whole period) ---\n');

[p_nppv_alfred, t_nppv_alf] = manual_ttest2(ts_nppv_alfred, ts_base_nppv);
[p_spm_alfred,  t_spm_alf]  = manual_ttest2(ts_spm_alfred,  ts_base_spm);
[p_nppv_fina,   t_nppv_fin] = manual_ttest2(ts_nppv_fina,   ts_base_nppv);
[p_spm_fina,    t_spm_fin]  = manual_ttest2(ts_spm_fina,    ts_base_spm);

[p_anom_nppv_alfred, t_anom_nppv_alf] = manual_ttest2(anom_ts_nppv_alfred, zeros(size(anom_ts_nppv_alfred)));
[p_anom_spm_alfred,  t_anom_spm_alf]  = manual_ttest2(anom_ts_spm_alfred,  zeros(size(anom_ts_spm_alfred)));
[p_anom_nppv_fina,   t_anom_nppv_fin] = manual_ttest2(anom_ts_nppv_fina,   zeros(size(anom_ts_nppv_fina)));
[p_anom_spm_fina,    t_anom_spm_fin]  = manual_ttest2(anom_ts_spm_fina,    zeros(size(anom_ts_spm_fina)));

disp(' ')
disp('============================================================')
disp('  T-TEST RESULTS — Time series in overlap box')
disp('============================================================')

% FIX (B): All four blocks now have complete fprintf output including SPM Fina
fprintf('\n--- ALFRED vs BASELINE ---\n')
fprintf('NPPV: t = %.3f, p = %.5f  -> %s\n', t_nppv_alf, p_nppv_alfred, interpret_p(p_nppv_alfred))
fprintf('SPM:  t = %.3f, p = %.5f  -> %s\n', t_spm_alf,  p_spm_alfred,  interpret_p(p_spm_alfred))

fprintf('\n--- FINA vs BASELINE ---\n')
fprintf('NPPV: t = %.3f, p = %.5f  -> %s\n', t_nppv_fin, p_nppv_fina, interpret_p(p_nppv_fina))
fprintf('SPM:  t = %.3f, p = %.5f  -> %s\n', t_spm_fin,  p_spm_fina,  interpret_p(p_spm_fina))

fprintf('\n--- ALFRED ANOMALY vs ZERO ---\n')
fprintf('NPPV anom: t = %.3f, p = %.5f  -> %s\n', t_anom_nppv_alf, p_anom_nppv_alfred, interpret_p(p_anom_nppv_alfred))
fprintf('SPM  anom: t = %.3f, p = %.5f  -> %s\n', t_anom_spm_alf,  p_anom_spm_alfred,  interpret_p(p_anom_spm_alfred))

fprintf('\n--- FINA ANOMALY vs ZERO ---\n')
fprintf('NPPV anom: t = %.3f, p = %.5f  -> %s\n', t_anom_nppv_fin, p_anom_nppv_fina, interpret_p(p_anom_nppv_fina))
fprintf('SPM  anom: t = %.3f, p = %.5f  -> %s\n', t_anom_spm_fin,  p_anom_spm_fina,  interpret_p(p_anom_spm_fina))

% Save t-test table to CSV
ttest_results = table( ...
    ["Alfred vs Baseline"; "Alfred vs Baseline"; "Fina vs Baseline"; "Fina vs Baseline"; ...
     "Alfred anomaly vs zero"; "Alfred anomaly vs zero"; "Fina anomaly vs zero"; "Fina anomaly vs zero"], ...
    ["NPPV"; "SPM"; "NPPV"; "SPM"; "NPPV anomaly"; "SPM anomaly"; "NPPV anomaly"; "SPM anomaly"], ...
    [t_nppv_alf; t_spm_alf; t_nppv_fin; t_spm_fin; ...
     t_anom_nppv_alf; t_anom_spm_alf; t_anom_nppv_fin; t_anom_spm_fin], ...
    [p_nppv_alfred; p_spm_alfred; p_nppv_fina; p_spm_fina; ...
     p_anom_nppv_alfred; p_anom_spm_alfred; p_anom_nppv_fina; p_anom_spm_fina], ...
    [p_nppv_alfred<0.05; p_spm_alfred<0.05; p_nppv_fina<0.05; p_spm_fina<0.05; ...
     p_anom_nppv_alfred<0.05; p_anom_spm_alfred<0.05; p_anom_nppv_fina<0.05; p_anom_spm_fina<0.05], ...
    'VariableNames', {'Test','Variable','t_value','p_value','Significant_005'});

csv_name = sprintf('ttest_results_NPPV_SPM_%s.csv', datestr(now,'yyyymmdd_HHMMSS'));
try
    writetable(ttest_results, fullfile(output_dir, csv_name));
    fprintf('\nT-test table saved as: %s\n', csv_name);
catch ME
    warning('Could not save CSV: %s\nClose the file in Excel and re-run if needed.', ME.message);
end

%% =========================================================================
%  5. FIGURES 
%
%  FIX (C): All time series use "days since start of cyclone period" on the
%  x-axis (day 0 = first day of data, e.g. 1 Feb for Alfred). The peak is
%  drawn as a vertical line at the correct offset from the start date.
%      days_alf = days(t_alf - t_alf(1))          % starts at 0
%      peak_day_alf = days(alfred_peak - alfred_start)   % peak line position
% =========================================================================

fprintf('\n--- Generating figures ---\n');

% Shared colour limits (Alfred & Fina on same scale for fair comparison)
clims.nppv_mean = safe_clim([0, max([alfred_nppv_mean(:); fina_nppv_mean(:)], [], 'omitnan')]);
clims.spm_mean  = safe_clim([0, max([alfred_spm_mean(:);  fina_spm_mean(:)],  [], 'omitnan')]);
clims.nppv_std  = safe_clim([0, max([alfred_nppv_std(:);  fina_nppv_std(:)],  [], 'omitnan')]);
clims.spm_std   = safe_clim([0, max([alfred_spm_std(:);   fina_spm_std(:)],   [], 'omitnan')]);

nppv_anom_abs = percentile_no_toolbox(abs([peak3_anom_nppv_alfred(:); peak3_anom_nppv_fina(:)]), 98);
spm_anom_abs  = percentile_no_toolbox(abs([peak3_anom_spm_alfred(:);  peak3_anom_spm_fina(:)]),  98);
clims.nppv_anom = safe_clim([-nppv_anom_abs  nppv_anom_abs]);
clims.spm_anom  = safe_clim([-spm_anom_abs   spm_anom_abs]);

% =====================================================================
% FIGURE 0 — NPPV & SPM: Baseline time series
% =====================================================================
plot_baseline_timeseries( ...
    all_time_nppv_base, ts_base_nppv, ...
    all_time_spm_base,  ts_base_spm, ...
    'NPPV and SPM during Baseline period (Cyclone-Free)', ...
    'Baseline_NPPV_SPM_timeseries', output_dir);
% =====================================================================
% FIGURE 1 — NPPV: Time series in overlap box (raw values)
% X-axis = days since start of cyclone period (day 0 = start, not peak)
% =====================================================================
plot_dual_cyclone_timeseries( ...
    t_nppv_alfred, ts_nppv_alfred, alfred_start, alfred_peak, ...
    t_nppv_fina,   ts_nppv_fina,   fina_start,   fina_peak, ...
    all_time_nppv_base, ts_base_nppv, ...
    'NPPV — Time Series in Overlap Box', ...
    'NPPV (mg C m^{-2} d^{-1})', ...
    'NPPV_timeseries_box', output_dir, false);

% =====================================================================
% FIGURE 2 — SPM: Time series in overlap box (raw values)
% =====================================================================
plot_dual_cyclone_timeseries( ...
    t_spm_alfred, ts_spm_alfred, alfred_start, alfred_peak, ...
    t_spm_fina,   ts_spm_fina,   fina_start,   fina_peak, ...
    all_time_spm_base, ts_base_spm, ...
    'SPM — Time Series in Overlap Box', ...
    'SPM (g m^{-3})', ...
    'SPM_timeseries_box', output_dir, false);

% =====================================================================
% FIGURE 3 — NPPV: Anomaly time series in overlap box
% =====================================================================
plot_dual_anomaly_timeseries( ...
    t_nppv_alfred, anom_ts_nppv_alfred, alfred_start, alfred_peak, ...
    t_nppv_fina,   anom_ts_nppv_fina,   fina_start,   fina_peak, ...
    'NPPV — Anomaly Time Series in Overlap Box', ...
    'ΔNPPV (mg C m^{-2} d^{-1})', ...
    'NPPV_anomaly_timeseries_box', output_dir);

% =====================================================================
% FIGURE 4 — SPM: Anomaly time series in overlap box
% =====================================================================
plot_dual_anomaly_timeseries( ...
    t_spm_alfred, anom_ts_spm_alfred, alfred_start, alfred_peak, ...
    t_spm_fina,   anom_ts_spm_fina,   fina_start,   fina_peak, ...
    'SPM — Anomaly Time Series in Overlap Box', ...
    'ΔSPM (g m^{-3})', ...
    'SPM_anomaly_timeseries_box', output_dir);

% =====================================================================
% FIGURE 5 — NPPV: Mean & Std maps (Alfred vs Fina)
% =====================================================================
plot_mean_std_comparison( ...
    alfred_nppv_mean, alfred_nppv_std, ...
    fina_nppv_mean,   fina_nppv_std, ...
    lat_nppv, lon_nppv, lat_box, lon_box, ...
    'NPPV Mean and Standard Deviation', ...
    'NPPV_Mean_STD_Comparison', ...
    'mg C m^{-2} d^{-1}', output_dir, ...
    clims.nppv_mean, clims.nppv_std, LAND_COLOR);

% =====================================================================
% FIGURE 6 — SPM: Mean & Std maps (Alfred vs Fina)
% =====================================================================
plot_mean_std_comparison( ...
    alfred_spm_mean, alfred_spm_std, ...
    fina_spm_mean,   fina_spm_std, ...
    lat_spm, lon_spm, lat_box, lon_box, ...
    'SPM Mean and Standard Deviation', ...
    'SPM_Mean_STD_Comparison', ...
    'g m^{-3}', output_dir, ...
    clims.spm_mean, clims.spm_std, LAND_COLOR);

% =====================================================================
% FIGURE 7 — NPPV: 3-day anomaly maps (total + significant)
% =====================================================================
plot_anomaly_comparison( ...
    peak3_anom_nppv_alfred, sig_anom_nppv_alfred, ...
    peak3_anom_nppv_fina,   sig_anom_nppv_fina, ...
    lat_nppv, lon_nppv, lat_box, lon_box, ...
    'NPPV 3-day Anomaly around Peak (vs Baseline)', ...
    'NPPV_3day_Anomaly_Comparison', ...
    'mg C m^{-2} d^{-1}', output_dir, clims.nppv_anom, LAND_COLOR, NONSIG_COLOR);

% =====================================================================
% FIGURE 8 — SPM: 3-day anomaly maps (total + significant)
% =====================================================================
plot_anomaly_comparison( ...
    peak3_anom_spm_alfred, sig_anom_spm_alfred, ...
    peak3_anom_spm_fina,   sig_anom_spm_fina, ...
    lat_spm, lon_spm, lat_box, lon_box, ...
    'SPM 3-day Anomaly around Peak (vs Baseline)', ...
    'SPM_3day_Anomaly_Comparison', ...
    'g m^{-3}', output_dir, clims.spm_anom, LAND_COLOR, NONSIG_COLOR);

%% =========================================================================
%  6. SAVE RESULTS
% =========================================================================

save(fullfile(output_dir, 'NPPV_SPM_results.mat'), ...
    'lon_nppv', 'lat_nppv', 'lon_spm', 'lat_spm', ...
    'base_nppv_mean', 'base_nppv_std', 'base_spm_mean', 'base_spm_std', ...
    'alfred_nppv_mean', 'alfred_nppv_std', 'alfred_spm_mean', 'alfred_spm_std', ...
    'fina_nppv_mean',   'fina_nppv_std',   'fina_spm_mean',   'fina_spm_std', ...
    'anom_nppv_alfred', 'anom_nppv_fina',  'anom_spm_alfred', 'anom_spm_fina', ...
    'zscore_nppv_alfred','zscore_nppv_fina','zscore_spm_alfred','zscore_spm_fina', ...
    'peak3_anom_nppv_alfred', 'sig_anom_nppv_alfred', 'pmap_nppv_alfred', ...
    'peak3_anom_spm_alfred',  'sig_anom_spm_alfred',  'pmap_spm_alfred', ...
    'peak3_anom_nppv_fina',   'sig_anom_nppv_fina',   'pmap_nppv_fina', ...
    'peak3_anom_spm_fina',    'sig_anom_spm_fina',    'pmap_spm_fina', ...
    'ts_nppv_alfred', 't_nppv_alfred', 'ts_spm_alfred', 't_spm_alfred', ...
    'ts_nppv_fina',   't_nppv_fina',   'ts_spm_fina',   't_spm_fina');

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Figures and results saved in: %s\n', output_dir);

%% =========================================================================
%  7. LOCAL FUNCTIONS
% =========================================================================

% -------------------------------------------------------------------------
function [data_cube, time_out, lon_out, lat_out] = load_baseline_data( ...
        files, var_name, lon_var, lat_var, time_var, is_nppv)
% Loads and concatenates baseline data, excluding forbidden cyclone periods.
% is_nppv = true  → data has 4 dims [lon x lat x depth x time], read depth=1
% is_nppv = false → data has 3 dims [lon x lat x time]
%
% FIX (E): exclude periods are applied here via get_exclude_periods()

    exclude = get_exclude_periods();
    data_cube = [];
    time_out  = datetime.empty(0,1);
    time_out.TimeZone = 'UTC';
    lon_out   = [];
    lat_out   = [];

    for i = 1:numel(files)
        f = files{i};

        % Read coordinates once
        if isempty(lon_out)
            lon_out = double(ncread(f, lon_var));
            lat_out = double(ncread(f, lat_var));
            if isrow(lon_out); lon_out = lon_out'; end
            if isrow(lat_out); lat_out = lat_out'; end
        end

        % Read time
        raw_t  = ncread(f, time_var);
        t_unit = ncreadatt(f, time_var, 'units');
        t_all  = convert_cmems_time(raw_t, t_unit);

        % Keep only Nov-Mar months (Australian summer baseline)
        keep = month(t_all) >= 11 | month(t_all) <= 3;

        % Remove all forbidden cyclone periods
        for j = 1:size(exclude, 1)
            keep = keep & ~(t_all >= exclude(j,1) & t_all <= exclude(j,2));
        end

        idx = find(keep);
        if isempty(idx); continue; end

        % Read data
        if is_nppv
            % [lon x lat x depth x time] → surface only (depth index 1)
            raw = double(ncread(f, var_name, ...
                [1, 1, 1, idx(1)], [Inf, Inf, 1, numel(idx)]));
            raw = squeeze(raw);   % → [lon x lat x time]
        else
            raw = double(ncread(f, var_name, ...
                [1, 1, idx(1)], [Inf, Inf, numel(idx)]));
        end

        % Replace _FillValue with NaN
        try
            fv = ncreadatt(f, var_name, '_FillValue');
            raw(raw == fv) = NaN;
        catch; end

        % Secondary sentinel cleanup
        sentinel_candidates = [100, 1000, 9999, 99999, -999, -9999];
        for sv = sentinel_candidates
            n_hits = sum(raw(:) == sv);
            if n_hits > 0 && n_hits < 0.01 * numel(raw)
                raw(raw == sv) = NaN;
            end
        end

        data_cube = cat(3, data_cube, raw);
        time_out  = [time_out; t_all(keep)]; %#ok<AGROW>
    end
end

% -------------------------------------------------------------------------
function [data_cube, time_out] = load_cyclone_data( ...
        file, var_name, time_var, t_start, t_end, is_nppv)
% Loads data for a single cyclone period.

    raw_t  = ncread(file, time_var);
    t_unit = ncreadatt(file, time_var, 'units');
    t_all  = convert_cmems_time(raw_t, t_unit);

    mask = (t_all >= t_start) & (t_all <= t_end);
    idx  = find(mask);
    if isempty(idx)
        error('No timesteps found in file %s for the given date range.', file);
    end

    if is_nppv
        raw = double(ncread(file, var_name, ...
            [1, 1, 1, idx(1)], [Inf, Inf, 1, numel(idx)]));
        raw = squeeze(raw);
    else
        raw = double(ncread(file, var_name, ...
            [1, 1, idx(1)], [Inf, Inf, numel(idx)]));
    end

    try
        fv = ncreadatt(file, var_name, '_FillValue');
        raw(raw == fv) = NaN;
    catch; end

    % Secondary fill-value check
    sentinel_candidates = [100, 1000, 9999, 99999, -999, -9999];
    for sv = sentinel_candidates
        n_hits = sum(raw(:) == sv);
        if n_hits > 0 && n_hits < 0.01 * numel(raw)
            fprintf('  [INFO] Masking %d pixel(s) with sentinel value %.0f in %s\n', ...
                n_hits, sv, var_name);
            raw(raw == sv) = NaN;
        end
    end

    data_cube = raw;
    time_out  = t_all(mask);
end

% -------------------------------------------------------------------------
function ts = extract_box_ts(data_cube, lon, lat, lon_box, lat_box)
% Computes spatial mean within a lon/lat box for each timestep.
% data_cube: [lon x lat x time]

    lon_i = lon >= lon_box(1) & lon <= lon_box(2);
    lat_i = lat >= lat_box(1) & lat <= lat_box(2);

    nt = size(data_cube, 3);
    ts = NaN(nt, 1);
    for k = 1:nt
        frame = data_cube(:,:,k);
        sub   = frame(lon_i, lat_i);
        ts(k) = mean(sub(:), 'omitnan');
    end
end

% -------------------------------------------------------------------------
function base_ts = interp_baseline_to_event(t_event, t_base, ts_base)
% Maps baseline time series to the event dates using season-day matching.
% Works across year boundaries (e.g. Dec 2011 – Jan 2012).

    base_ts = NaN(numel(t_event), 1);

    ref_nov1 = datetime(2001, 11, 1, 'TimeZone', 'UTC');

    % Compute season-day index for each baseline date
    season_day_base = NaN(numel(t_base), 1);
    for k = 1:numel(t_base)
        m = month(t_base(k)); d = day(t_base(k));
        if m >= 11
            season_day_base(k) = days(datetime(2001,m,d,'TimeZone','UTC') - ref_nov1) + 1;
        else
            season_day_base(k) = days(datetime(2002,m,d,'TimeZone','UTC') - ref_nov1) + 1;
        end
    end

    for i = 1:numel(t_event)
        m = month(t_event(i)); d = day(t_event(i));
        if m >= 11
            sd = days(datetime(2001,m,d,'TimeZone','UTC') - ref_nov1) + 1;
        else
            sd = days(datetime(2002,m,d,'TimeZone','UTC') - ref_nov1) + 1;
        end

        % Average all baseline values on that season day
        match = season_day_base == sd;
        if any(match)
            base_ts(i) = mean(ts_base(match), 'omitnan');
        end
    end
end

% -------------------------------------------------------------------------
function peak3 = build_3day_peak_map(data_cube, time_vec, peak_date)
% Returns the mean map over [peak-1, peak, peak+1] days.

    idx = isbetween(time_vec, peak_date - days(1), peak_date + days(1));
    if ~any(idx)
        error('No data found within 3 days of peak date %s.', datestr(peak_date));
    end
    peak3 = mean(data_cube(:,:,idx), 3, 'omitnan');
end

% -------------------------------------------------------------------------
function stack = build_baseline_3day_stack(data_cube, time_vec, event_peak)
% Builds a stack of baseline maps for the same 3 calendar days as event peak,
% across all baseline years. Each year contributes one mean map (layer).

    exclude = get_exclude_periods();
    target_dates = event_peak + days(-1:1);   % [peak-1, peak, peak+1]
    stack = [];

    years_in_base = unique(year(time_vec));

    for yr = years_in_base'
        keep = false(size(time_vec));
        for d = 1:numel(target_dates)
            keep = keep | ( month(time_vec) == month(target_dates(d)) & ...
                            day(time_vec)   == day(target_dates(d))   & ...
                            year(time_vec)  == yr );
        end

        % Remove forbidden periods
        for j = 1:size(exclude, 1)
            keep = keep & ~(time_vec >= exclude(j,1) & time_vec <= exclude(j,2));
        end

        if any(keep)
            yr_mean = mean(data_cube(:,:,keep), 3, 'omitnan');
            stack   = cat(3, stack, yr_mean);
        end
    end

    if isempty(stack)
        warning('baseline_3day_stack: no matching days found for peak %s', datestr(event_peak));
        sz    = size(data_cube);
        stack = NaN(sz(1), sz(2), 1);
    end
end

% -------------------------------------------------------------------------
function [sig_map, pmap] = gridcell_ttest_mask(cyclone_peak3, base_stack, alpha)
% One-sample t-test at each grid cell:
%   H0: cyclone 3-day mean equals baseline mean
%   base_stack: [lon x lat x n_years]
%   cyclone_peak3: [lon x lat]

    base_mean = mean(base_stack, 3, 'omitnan');
    base_std  = std(base_stack,  0, 3, 'omitnan');
    n         = sum(~isnan(base_stack), 3);

    se   = base_std ./ sqrt(n);
    tmap = (cyclone_peak3 - base_mean) ./ se;

    % Two-tailed p-value (no Statistics Toolbox needed)
    pmap = erfc(abs(tmap) ./ sqrt(2));

    sig_mask = pmap < alpha;
    sig_mask(isnan(pmap)) = false;
    sig_mask(n < 2)       = false;
    sig_mask(se == 0)     = false;

    sig_map = cyclone_peak3;
    sig_map(~sig_mask) = NaN;
end

% -------------------------------------------------------------------------
function [p, t] = manual_ttest2(x, y)
% Two-sample Welch t-test, no Statistics Toolbox needed.

    x = x(~isnan(x));  y = y(~isnan(y));
    nx = numel(x);     ny = numel(y);
    if nx < 2 || ny < 2;  p = NaN; t = NaN;  return;  end

    mx = mean(x);  my = mean(y);
    vx = var(x);   vy = var(y);
    se = sqrt(vx/nx + vy/ny);

    if se == 0 || isnan(se);  p = NaN; t = NaN;  return;  end

    t = (mx - my) / se;
    p = erfc(abs(t) / sqrt(2));
end

% -------------------------------------------------------------------------
function t = convert_cmems_time(raw_time, units_str)
% Converts CMEMS numeric time to MATLAB datetime.

    raw_time  = double(raw_time);
    units_str = strtrim(units_str);

    if contains(units_str, 'seconds since 1970-01-01')
        t = datetime(raw_time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    elseif contains(units_str, 'hours since')
        ref = datetime(strrep(units_str,'hours since ',''), ...
                       'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','UTC');
        t = ref + hours(raw_time);
    elseif contains(units_str, 'days since')
        ref = datetime(strrep(units_str,'days since ',''), ...
                       'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','UTC');
        t = ref + days(raw_time);
    elseif contains(units_str, 'seconds since')
        ref = datetime(strrep(units_str,'seconds since ',''), ...
                       'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','UTC');
        t = ref + seconds(raw_time);
    else
        error('Unknown time unit: %s', units_str);
    end
    t.TimeZone = 'UTC';
end

% -------------------------------------------------------------------------
function exclude = get_exclude_periods()
% All cyclone-contaminated periods to exclude from baseline.
% 8 periods total.
    exclude = [
        datetime(2022, 2,  7,'TimeZone','UTC')  datetime(2022, 2,12,'TimeZone','UTC');  % Dovi
        datetime(2021,12,  9,'TimeZone','UTC')  datetime(2021,12,15,'TimeZone','UTC');  % Ruby
        datetime(2021, 2, 27,'TimeZone','UTC')  datetime(2021, 3, 5,'TimeZone','UTC');  % Niran
        datetime(2021, 1, 24,'TimeZone','UTC')  datetime(2021, 2, 1,'TimeZone','UTC');  % Lucas
        datetime(2021, 1, 16,'TimeZone','UTC')  datetime(2021, 1,19,'TimeZone','UTC');  % Kimi
        datetime(2020, 4,  1,'TimeZone','UTC')  datetime(2020, 4,10,'TimeZone','UTC');  % Harold
        datetime(2020, 3, 14,'TimeZone','UTC')  datetime(2020, 3,15,'TimeZone','UTC');  % Gretel
        datetime(2020, 2,  4,'TimeZone','UTC')  datetime(2020, 2,14,'TimeZone','UTC');  % Uesi
    ];
end

% -------------------------------------------------------------------------
%  FIGURE FUNCTIONS  
% -------------------------------------------------------------------------

function y_out = fill_gaps(x, y)
% Fills NaN gaps in y using linear interpolation over x (for plotting only).
    y = y(:); x = x(:);
    valid = ~isnan(y);
    if sum(valid) < 2
        y_out = y;
        return
    end
    y_out = interp1(x(valid), y(valid), x, 'linear', NaN);
end

% -------------------------------------------------------------------------
function plot_dual_cyclone_timeseries( ...
        t_alf, ts_alf, start_alf, peak_alf, ...
        t_fin, ts_fin, start_fin, peak_fin, ...
        t_base, ts_base, ...
        title_str, ylabel_str, filename, output_dir, show_zero)
%
% X-axis = days since start of cyclone period (day 0 = first data point).
% Peak is drawn as a vertical line at its correct offset from the start.

    figure('Position', [100 100 900 500]);

    % Days since start of each cyclone period (both start at 0)
    days_alf = days(t_alf - start_alf);
    days_fin = days(t_fin - start_fin);

    % Position of peak line relative to start
    peak_day_alf = days(peak_alf - start_alf);
    peak_day_fin = days(peak_fin - start_fin);

    % Fill NaN gaps with linear interpolation (display only)
    ts_alf_plot = fill_gaps(days_alf, ts_alf);
    ts_fin_plot = fill_gaps(days_fin, ts_fin);

    base_mean_val = mean(ts_base, 'omitnan');
    base_std_val  = std(ts_base,  'omitnan');

    ax = gca;
    hold(ax, 'on');

    h1 = plot(ax, days_alf, ts_alf_plot, '-', ...
        'Color', [0.20 0.50 0.90], 'LineWidth', 2.0, 'DisplayName', 'Cyclone Alfred');
    h2 = plot(ax, days_fin, ts_fin_plot, '-', ...
        'Color', [0.13 0.67 0.42], 'LineWidth', 2.0, 'DisplayName', 'Cyclone Fina');

    yline(ax, base_mean_val, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
        'Label', sprintf('Baseline mean = %.3f', base_mean_val), ...
        'LabelHorizontalAlignment', 'left');
    yline(ax, base_mean_val + base_std_val, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
    yline(ax, base_mean_val - base_std_val, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);

    if show_zero
        yline(ax, 0, '--r', 'LineWidth', 1.2);
    end

    % Peak lines at correct day offset from start (one per cyclone)
    xline(ax, peak_day_alf, '--', 'Color', [0.20 0.50 0.90], 'LineWidth', 1.8, ...
        'Label', sprintf('Alfred peak (day %d)', round(peak_day_alf)), ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
    xline(ax, peak_day_fin, '--', 'Color', [0.13 0.67 0.42], 'LineWidth', 1.8, ...
        'Label', sprintf('Fina peak (day %d)', round(peak_day_fin)), ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

    legend(ax, [h1 h2], {'Cyclone Alfred','Cyclone Fina'}, ...
        'Location', 'best', 'FontSize', 9);
    xlabel(ax, 'Days since start of cyclone period (day 0 = start)', 'FontSize', 11);
    ylabel(ax, ylabel_str, 'FontSize', 11, 'FontWeight', 'bold');
    title(ax, title_str, 'FontSize', 13, 'FontWeight', 'bold');
    set(ax, 'FontSize', 10, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
        'GridAlpha', 0.25, 'GridLineStyle', '--');
    hold(ax, 'off');

    exportgraphics(gcf, fullfile(output_dir, [filename '.png']), 'Resolution', 300);
    exportgraphics(gcf, fullfile(output_dir, [filename '.pdf']), 'ContentType', 'vector');
    close(gcf);
end

% -------------------------------------------------------------------------
function plot_dual_anomaly_timeseries( ...
        t_alf, anom_alf, start_alf, peak_alf, ...
        t_fin, anom_fin, start_fin, peak_fin, ...
        title_str, ylabel_str, filename, output_dir)
%
% X-axis = days since start of cyclone period (day 0 = first data point).
% Peak is drawn as a vertical line at its correct offset from the start.

    figure('Position', [100 100 900 500]);

    % Days since start of each cyclone period (both start at 0)
    days_alf = days(t_alf - start_alf);
    days_fin = days(t_fin - start_fin);

    % Position of peak line relative to start
    peak_day_alf = days(peak_alf - start_alf);
    peak_day_fin = days(peak_fin - start_fin);

    anom_alf_plot = fill_gaps(days_alf, anom_alf);
    anom_fin_plot = fill_gaps(days_fin, anom_fin);

    ax = gca;
    hold(ax, 'on');

    h1 = plot(ax, days_alf, anom_alf_plot, '-', ...
        'Color', [0.20 0.50 0.90], 'LineWidth', 2.0, 'DisplayName', 'Alfred anomaly');
    h2 = plot(ax, days_fin, anom_fin_plot, '-', ...
        'Color', [0.13 0.67 0.42], 'LineWidth', 2.0, 'DisplayName', 'Fina anomaly');

    yline(ax, 0, '--k', 'LineWidth', 1.2);

    % Peak lines at correct day offset from start (one per cyclone)
    xline(ax, peak_day_alf, '--', 'Color', [0.20 0.50 0.90], 'LineWidth', 1.8, ...
        'Label', sprintf('Alfred peak (day %d)', round(peak_day_alf)), ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
    xline(ax, peak_day_fin, '--', 'Color', [0.13 0.67 0.42], 'LineWidth', 1.8, ...
        'Label', sprintf('Fina peak (day %d)', round(peak_day_fin)), ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

    legend(ax, [h1 h2], {'Alfred anomaly','Fina anomaly'}, ...
        'Location', 'best', 'FontSize', 9);
    xlabel(ax, 'Days since start of cyclone period (day 0 = start)', 'FontSize', 11);
    ylabel(ax, ylabel_str, 'FontSize', 11, 'FontWeight', 'bold');
    title(ax, title_str, 'FontSize', 13, 'FontWeight', 'bold');
    set(ax, 'FontSize', 10, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
        'GridAlpha', 0.25, 'GridLineStyle', '--');
    hold(ax, 'off');

    exportgraphics(gcf, fullfile(output_dir, [filename '.png']), 'Resolution', 300);
    exportgraphics(gcf, fullfile(output_dir, [filename '.pdf']), 'ContentType', 'vector');
    close(gcf);
end

% -------------------------------------------------------------------------
function plot_mean_std_comparison( ...
        mean_a, std_a, mean_f, std_f, ...
        lat, lon, lat_box, lon_box, ...
        title_str, filename, var_label, output_dir, ...
        mean_lim, std_lim, land_color)

    figure('Position', [100 100 1100 800]);
    cmap = turbo(256);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile
    plot_single_map(mean_a, lat, lon, lat_box, lon_box, cmap, land_color, mean_lim);
    title('Alfred — Mean');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    nexttile
    plot_single_map(std_a, lat, lon, lat_box, lon_box, cmap, land_color, std_lim);
    title('Alfred — Std Dev');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    nexttile
    plot_single_map(mean_f, lat, lon, lat_box, lon_box, cmap, land_color, mean_lim);
    title('Fina — Mean');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    nexttile
    plot_single_map(std_f, lat, lon, lat_box, lon_box, cmap, land_color, std_lim);
    title('Fina — Std Dev');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    sgtitle(title_str, 'FontSize', 14, 'FontWeight', 'bold');

    exportgraphics(gcf, fullfile(output_dir, [filename '.png']), 'Resolution', 300);
    exportgraphics(gcf, fullfile(output_dir, [filename '.pdf']), 'ContentType', 'vector');
    close(gcf);
end

% -------------------------------------------------------------------------
function plot_anomaly_comparison( ...
        anom_a, sig_a, anom_f, sig_f, ...
        lat, lon, lat_box, lon_box, ...
        title_str, filename, var_label, output_dir, ...
        anom_lim, land_color, nonsig_color)

    figure('Position', [100 100 1100 800]);
    cmap = redblue(256);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile
    plot_single_map(anom_a, lat, lon, lat_box, lon_box, cmap, land_color, anom_lim);
    title('Alfred — 3-day mean anomaly');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    nexttile
    plot_single_map(sig_a, lat, lon, lat_box, lon_box, cmap, nonsig_color, anom_lim);
    title('Alfred — significant anomaly (p < 0.05)');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    nexttile
    plot_single_map(anom_f, lat, lon, lat_box, lon_box, cmap, land_color, anom_lim);
    title('Fina — 3-day mean anomaly');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    nexttile
    plot_single_map(sig_f, lat, lon, lat_box, lon_box, cmap, nonsig_color, anom_lim);
    title('Fina — significant anomaly (p < 0.05)');
    cb = colorbar; cb.Label.String = var_label; cb.Label.FontWeight = 'bold';

    sgtitle(title_str, 'FontSize', 14, 'FontWeight', 'bold');

    exportgraphics(gcf, fullfile(output_dir, [filename '.png']), 'Resolution', 300);
    exportgraphics(gcf, fullfile(output_dir, [filename '.pdf']), 'ContentType', 'vector');
    close(gcf);
end

% -------------------------------------------------------------------------
function plot_single_map(map_data, lat, lon, lat_box, lon_box, map_cmap, background_color, clim_values)
% Plots a 2D map using imagesc

    h = imagesc(lon, lat, map_data');
    set(gca, 'YDir', 'normal');
    set(gca, 'Color', background_color);
    set(h, 'AlphaData', ~isnan(map_data'));
    colormap(gca, map_cmap);
    caxis(clim_values);
    xlabel('Longitude (°E)', 'FontSize', 10);
    ylabel('Latitude (°N)',  'FontSize', 10);
    set(gca, 'FontSize', 10, 'Box', 'on');
    grid on; hold on;

    % Draw overlap box in red
    x_box = [lon_box(1) lon_box(2) lon_box(2) lon_box(1) lon_box(1)];
    y_box = [lat_box(1) lat_box(1) lat_box(2) lat_box(2) lat_box(1)];
    plot(x_box, y_box, 'r-', 'LineWidth', 2);
    hold off;
end

% -------------------------------------------------------------------------
function clim_out = safe_clim(clim_in)
    clim_out = clim_in;
    if any(isnan(clim_out)) || any(isinf(clim_out)); clim_out = [0 1]; end
    if clim_out(1) == clim_out(2); clim_out = clim_out + [-1 1]*eps; end
end

% -------------------------------------------------------------------------
function p = percentile_no_toolbox(x, perc)
    x = sort(x(~isnan(x)));
    if isempty(x); p = NaN; return; end
    idx = max(1, min(round((perc/100)*numel(x)), numel(x)));
    p = x(idx);
end

% -------------------------------------------------------------------------
function txt = interpret_p(p)
    if isnan(p);         txt = 'not available';
    elseif p < 0.01;     txt = 'highly significant (p < 0.01)';
    elseif p < 0.05;     txt = 'significant (p < 0.05)';
    else;                txt = 'not significant (p >= 0.05)';
    end
end

% -------------------------------------------------------------------------
function cmap = redblue(n)
    if nargin < 1; n = 256; end
    n1 = floor(n/2); n2 = n - n1;
    bot = [0 0 0.5]; mid = [1 1 1]; top = [0.5 0 0];
    cmap = [
        linspace(bot(1),mid(1),n1)' linspace(bot(2),mid(2),n1)' linspace(bot(3),mid(3),n1)';
        linspace(mid(1),top(1),n2)' linspace(mid(2),top(2),n2)' linspace(mid(3),top(3),n2)'
    ];
end
function plot_baseline_timeseries(t_nppv, ts_nppv, t_spm, ts_spm, ...
        title_str, filename, output_dir)

    figure('Position', [100 100 900 500]);
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ---- Top panel: NPPV ----
    ax1 = nexttile;
    hold(ax1, 'on');
    plot(ax1, t_nppv, ts_nppv, '-', 'Color', [0.20 0.50 0.90], ...
        'LineWidth', 1.5);
    yline(ax1, mean(ts_nppv,'omitnan'), '--', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.2, ...
        'Label', sprintf('Mean = %.3f', mean(ts_nppv,'omitnan')), ...
        'LabelHorizontalAlignment', 'left');
    ylabel(ax1, 'NPPV (mg C m^{-2} d^{-1})', 'FontSize', 10, ...
        'FontWeight', 'bold');
    set(ax1, 'FontSize', 10, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
        'GridAlpha', 0.25, 'GridLineStyle', '--');
    xtickformat(ax1, 'dd-MMM-yyyy');
    xtickangle(ax1, 30);
    hold(ax1, 'off');

    % ---- Bottom panel: SPM ----
    ax2 = nexttile;
    hold(ax2, 'on');
    plot(ax2, t_spm, ts_spm, '-', 'Color', [0.85 0.33 0.10], ...
        'LineWidth', 1.5);
    yline(ax2, mean(ts_spm,'omitnan'), '--', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.2, ...
        'Label', sprintf('Mean = %.3f', mean(ts_spm,'omitnan')), ...
        'LabelHorizontalAlignment', 'left');
    ylabel(ax2, 'SPM (g m^{-3})', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel(ax2, 'Date', 'FontSize', 10);
    set(ax2, 'FontSize', 10, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
        'GridAlpha', 0.25, 'GridLineStyle', '--');
    xtickformat(ax2, 'dd-MMM-yyyy');
    xtickangle(ax2, 30);
    hold(ax2, 'off');

    sgtitle(title_str, 'FontSize', 13, 'FontWeight', 'bold');

    exportgraphics(gcf, fullfile(output_dir, [filename '.png']), ...
        'Resolution', 300);
    exportgraphics(gcf, fullfile(output_dir, [filename '.pdf']), ...
        'ContentType', 'vector');
    close(gcf);
end
