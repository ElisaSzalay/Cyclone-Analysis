%% ===============================================================
% FINAL_ALL_IN_ONE_CHL.m
% FINAL VERSION (Baseline Time Series ONLY)
%% ===============================================================

clear; clc; close all;

addpath 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred\m_map'
cd 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred'

outputFolder='plots';
if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

%% FILES
file_raw='baseline_CH.nc';
file_out='baseline_CH_final.nc';
file_alf='Alfred_CH.nc';
file_fin='Fina_CH.nc';

%% READ BASELINE RAW
lon=double(ncread(file_raw,'longitude'));
lat=double(ncread(file_raw,'latitude'));

time=double(ncread(file_raw,'time'));
time=datetime(time,'ConvertFrom','posixtime');

CHL=double(ncread(file_raw,'CHL'));
CHL=squeeze(CHL);

%% FILTER DATES
m=month(time);
idx=(m==11)|(m==12)|(m==1)|(m==2)|(m==3);

ex=false(size(time));

ex=ex|(time>=datetime(2020,2,4)  & time<=datetime(2020,2,14,23,59,59));
ex=ex|(time>=datetime(2020,3,14) & time<=datetime(2020,3,15,23,59,59));
ex=ex|(time>=datetime(2021,1,16) & time<=datetime(2021,1,19,23,59,59));
ex=ex|(time>=datetime(2021,2,27) & time<=datetime(2021,3,5,23,59,59));
ex=ex|(time>=datetime(2021,12,9) & time<=datetime(2021,12,15,23,59,59));
ex=ex|(time>=datetime(2022,2,7) & time<=datetime(2022,2,12,23,59,59));

idx=idx & ~ex;

time_b=time(idx);
CHL_b=CHL(:,:,idx);

%% SAVE FILTERED BASELINE
if isfile(file_out); delete(file_out); end

nccreate(file_out,'longitude','Dimensions',{'longitude',length(lon)});
nccreate(file_out,'latitude','Dimensions',{'latitude',length(lat)});
nccreate(file_out,'time','Dimensions',{'time',length(time_b)});
nccreate(file_out,'CHL',...
'Dimensions',{'longitude',length(lon),'latitude',length(lat),'time',length(time_b)});

ncwrite(file_out,'longitude',lon);
ncwrite(file_out,'latitude',lat);
ncwrite(file_out,'time',posixtime(time_b));
ncwrite(file_out,'CHL',single(CHL_b));

%% READ CYCLONES
time_a=datetime(double(ncread(file_alf,'time')),'ConvertFrom','posixtime');
time_f=datetime(double(ncread(file_fin,'time')),'ConvertFrom','posixtime');

CHL_a=squeeze(double(ncread(file_alf,'CHL')));
CHL_f=squeeze(double(ncread(file_fin,'CHL')));

%% BOX
ilon=find(lon>=155 & lon<=156);
ilat=find(lat>=-19.5 & lat<=-18.5);

CHL_b=CHL_b(ilon,ilat,:);
CHL_a=CHL_a(ilon,ilat,:);
CHL_f=CHL_f(ilon,ilat,:);

lon_sub=lon(ilon);
lat_sub=lat(ilat);
[Lon,Lat]=meshgrid(lon_sub,lat_sub);

%% MAPS
mean_b=mean(CHL_b,3,'omitnan');
mean_a=mean(CHL_a,3,'omitnan');
mean_f=mean(CHL_f,3,'omitnan');

anom_a=mean_a-mean_b;
anom_f=mean_f-mean_b;

figure
m_proj('mercator','lon',[155 156],'lat',[-19.5 -18.5])
m_pcolor(Lon,Lat,mean_b'); shading interp
m_coast('patch',[0.7 0.7 0.7]); m_grid
colorbar; title('Baseline Mean CHL')
exportgraphics(gcf,fullfile(outputFolder,'CHL_BaselineMean.png'));

figure
m_proj('mercator','lon',[155 156],'lat',[-19.5 -18.5])
m_pcolor(Lon,Lat,anom_a'); shading interp
m_coast('patch',[0.7 0.7 0.7]); m_grid
colorbar; title('Alfred CHL Anomaly')
exportgraphics(gcf,fullfile(outputFolder,'CHL_AlfredAnomaly.png'));

figure
m_proj('mercator','lon',[155 156],'lat',[-19.5 -18.5])
m_pcolor(Lon,Lat,anom_f'); shading interp
m_coast('patch',[0.7 0.7 0.7]); m_grid
colorbar; title('Fina CHL Anomaly')
exportgraphics(gcf,fullfile(outputFolder,'CHL_FinaAnomaly.png'));

%% BASELINE TIME SERIES ONLY
ts_b=squeeze(mean(mean(CHL_b,1,'omitnan'),2,'omitnan'));

figure
plot(time_b,ts_b,'k','LineWidth',2)
grid on
title('Baseline Chlorophyll-a Time Series')
xlabel('Time')
ylabel('CHL')
exportgraphics(gcf,fullfile(outputFolder,'CHL_Baseline_TimeSeries.png'));

%% REGRESSION
ts_a=squeeze(mean(mean(CHL_a,1,'omitnan'),2,'omitnan'));
ts_f=squeeze(mean(mean(CHL_f,1,'omitnan'),2,'omitnan'));

n=min(length(ts_b),length(ts_a));
x=ts_b(1:n); y=ts_a(1:n);
p=polyfit(x,y,1); yfit=polyval(p,x); R=corrcoef(x,y); R2=R(1,2)^2;

figure
scatter(x,y,60,'filled'); hold on
plot(x,yfit,'r','LineWidth',2)
grid on
title(['Alfred vs Baseline CHL R^2=',num2str(R2,'%.3f')])
exportgraphics(gcf,fullfile(outputFolder,'CHL_Reg_Alfred.png'));

n=min(length(ts_b),length(ts_f));
x=ts_b(1:n); y=ts_f(1:n);
p=polyfit(x,y,1); yfit=polyval(p,x); R=corrcoef(x,y); R2=R(1,2)^2;

figure
scatter(x,y,60,'filled'); hold on
plot(x,yfit,'r','LineWidth',2)
grid on
title(['Fina vs Baseline CHL R^2=',num2str(R2,'%.3f')])
exportgraphics(gcf,fullfile(outputFolder,'CHL_Reg_Fina.png'));