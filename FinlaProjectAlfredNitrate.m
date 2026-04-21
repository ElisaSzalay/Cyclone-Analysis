%% Final Project: Alfred
% Install the toolboxes

addpath 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred\m_map'
addpath 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred\seawater_ver3_3.1'
addpath 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred\unixtime2mat'

%% Define path (folder) where the NetCDF files are stored

cd 'C:\Users\greta\Desktop\1S Subjects\Digital Ocean\FinalProjectAlfred'
outputFolder = 'plots';
%% Specify the NetCDF file to read
filename = 'AlfredNo3.nc';

%% Display content of the files --> check description of each variable

ncdisp(filename) % ncdisp shows dimensions, variables, and attributes of NetCDF

%% Read in the time and convert it

time = double(ncread(filename,'time'));
time = datetime(time,'ConvertFrom','posixtime');

timee = datevec(time); % if you want to convert to datevec, breaks time into year, month, day

% Create string dates for naming/saving files accordingly at end
[r c] = size(timee);

for i=1:r
    year(i,:) = num2str(timee(i,1));
    if timee(i,2)<=9
        month(i,:) = ['0' num2str(timee(i,2))];
    else
        month(i,:) = num2str(timee(i,2));
    end
    
    if timee(i,3)<=9
        day(i,:) = ['0' num2str(timee(i,3))];
    else
        day(i,:) = num2str(timee(i,3));
    end
end

clear r c timee i 

day = string(day);
month = string(month);
year = string(year);

%% Read in the variables (not in a loop, since we only have one NetCDF file)
% index (1,1) because there is only one file, otherhwise we would have to create a loop

lat = double(ncread(filename,'latitude'));  % latitude (°N)
lon = double(ncread(filename,'longitude'));  % longitude (°E)
no3= double(ncread(filename,'no3'));  % mass concentration of chlorophyll a in seawater

%% ------------------------------------------------------------
%% SURFACE LAYER ONLY (first depth level)
%% ------------------------------------------------------------

no3_surface = squeeze(no3(:,:,1,:));

disp('--- SIZE SURFACE ---')
size(no3_surface)


%% ------------------------------------------------------------
%% MEAN AND STD THROUGH TIME
%% ------------------------------------------------------------

no3_mean = mean(no3_surface,3,'omitnan');
no3_std  = std(no3_surface,0,3,'omitnan');

disp('--- SIZE MEAN ---')
size(no3_mean)

disp('--- SIZE STD ---')
size(no3_std)

% Force 2D matrices
no3_mean = reshape(no3_mean,length(lon),length(lat));
no3_std  = reshape(no3_std ,length(lon),length(lat));

%% ------------------------------------------------------------
%% GRID
%% ------------------------------------------------------------

[Lon,Lat] = meshgrid(lon,lat);

%% ------------------------------------------------------------
%% MEAN MAP
%% ------------------------------------------------------------

figure

m_proj('mercator',...
       'lon',[min(lon) max(lon)],...
       'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,no3_mean')
shading interp

m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12,'tickdir','out')

colorbar
colormap(jet)

title('Mean Surface Nitrate - Cyclone Alfred')
exportgraphics(gcf, fullfile(outputFolder, 'Mean_Surface_Nitrate_Alfred.png'));

%% ------------------------------------------------------------
%% STD MAP
%% ------------------------------------------------------------

figure

m_proj('mercator',...
       'lon',[min(lon) max(lon)],...
       'lat',[min(lat) max(lat)])

m_pcolor(Lon,Lat,no3_std')
shading interp

m_coast('patch',[0.7 0.7 0.7],'edgecolor','k')
m_grid('fontsize',12,'tickdir','out')

colorbar
colormap(parula)

title('STD Surface Nitrate - Cyclone Alfred')
exportgraphics(gcf, fullfile(outputFolder, 'STD_Surface_Nitrate_Alfred.png'));