%% Final Project: Alfred
% WORKING VERSION
% Surface Nitrate (NO3)
% Analysis only for box:
% Longitude = 155 to 156 E
% Latitude  = -19.5 to -18.5 S

%% =====================================================
% TOOLBOXES
%% =====================================================

addpath 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred\m_map'
addpath 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred\seawater_ver3_3.1'
addpath 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred\unixtime2mat'

%% =====================================================
% PATH
%% =====================================================

cd 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred'
outputFolder = 'plots';

%% =====================================================
% FILE
%% =====================================================

filename = 'Alfred_NO3.nc';

%% =====================================================
% READ DATA
%% =====================================================

ncdisp(filename)

time = double(ncread(filename,'time'));
time = datetime(time,'ConvertFrom','posixtime');

lat = double(ncread(filename,'latitude'));
lon = double(ncread(filename,'longitude'));
no3 = double(ncread(filename,'no3'));

%% =====================================================
% SURFACE LAYER ONLY
%% =====================================================

NO3 = squeeze(no3(:,:,1,:));   % surface depth only

%% =====================================================
% SELECT TARGET BOX
%% =====================================================

ilon = find(lon >= 155 & lon <= 156);
ilat = find(lat >= -19.5 & lat <= -18.5);

lon_sub = lon(ilon);
lat_sub = lat(ilat);

NO3_sub = NO3(ilon,ilat,:);

%% =====================================================
% GRID
%% =====================================================

[Lon,Lat] = meshgrid(lon_sub,lat_sub);

%% =====================================================
% MEAN MAP
%% =====================================================

NO3_mean = mean(NO3_sub,3,'omitnan');

figure

m_proj('mercator',...
       'lon',[155 156],...
       'lat',[-19.5 -18.5])

m_pcolor(Lon,Lat,NO3_mean')
shading interp
hold on

box_lon = [155 156 156 155 155];
box_lat = [-18.5 -18.5 -19.5 -19.5 -18.5];

m_line(box_lon,box_lat,'color','r','LineWidth',2)

m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12,'tickdir','out')

colorbar
colormap(jet)

title('Mean Surface Nitrate (Alfred Box)')

exportgraphics(gcf,fullfile(outputFolder,'BOX_Mean_NO3_Alfred.png'));

%% =====================================================
% STD MAP
%% =====================================================

NO3_std = std(NO3_sub,0,3,'omitnan');

figure

m_proj('mercator',...
       'lon',[155 156],...
       'lat',[-19.5 -18.5])

m_pcolor(Lon,Lat,NO3_std')
shading interp
hold on

m_line(box_lon,box_lat,'color','r','LineWidth',2)

m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12,'tickdir','out')

colorbar
colormap(parula)

title('STD Surface Nitrate (Alfred Box)')

exportgraphics(gcf,fullfile(outputFolder,'BOX_STD_NO3_Alfred.png'));

%% =====================================================
% AREA-AVERAGED TIME SERIES
%% =====================================================

nt = length(time);
ts = nan(nt,1);

for t = 1:nt
    temp = NO3_sub(:,:,t);
    ts(t) = mean(temp(:),'omitnan');
end

%% =====================================================
% ANOMALY
%% =====================================================

ts_mean = mean(ts,'omitnan');
anomaly = ts - ts_mean;

[~,idx] = max(abs(anomaly));

fprintf('\n===== STRONGEST BOX ANOMALY (ALFRED NO3) =====\n');
fprintf('Value: %.4f\n', anomaly(idx));
fprintf('Date : %s\n', string(time(idx)));

%% =====================================================
% TIME SERIES
%% =====================================================

figure
plot(time,ts,'b','LineWidth',1.8)
grid on

title('Box-Averaged Surface Nitrate Time Series (Alfred)')
xlabel('Time')
ylabel('NO3')

exportgraphics(gcf,fullfile(outputFolder,'BOX_TimeSeries_NO3_Alfred.png'));

%% =====================================================
% ANOMALY PLOT
%% =====================================================

figure
plot(time,anomaly,'k','LineWidth',1.8)
hold on

yline(0,'r--')
scatter(time(idx),anomaly(idx),80,'r','filled')

grid on

title('Box-Averaged Surface Nitrate Anomaly (Alfred)')
xlabel('Time')
ylabel('Anomaly')

exportgraphics(gcf,fullfile(outputFolder,'BOX_Anomaly_NO3_Alfred.png'));