function [Year, Month, Day] = Data_Format_Interpolation_fast_pa(file, base, ext)
%{
If formated as an .xlsx and or .csv file, the input data should be compiled
where each column represents:
Time | Depth | Latitude | Longitude | Temperature | Ax | Ay | Az |
If formatted as a mat file, each variable named as follows:
time | depth |   lat    |   lon     |    temp     | ax | ay | az |
If a user does not have a certain variable the respective column should be
left blank, but the order/spacing should be maintained.
%}

%% Initialization
%month/day characteristics - reusable metrics
monthStartDay = [1 32 60 91 121 152 182 213 244 274 305 335;... %start day of each month
    1 32 61 92 122 153 183 214 245 275 306 336]; %row 1: Common Year row 2: Leap Year

%% Read In
switch ext %load in data (can be in .mat, .xlsx, or .csv)
    case 'mat'; load(file);
    case 'xlsx'; Data = xlsread(file);
    case 'csv'; Data = csvread(file);
end

%ensure variables possess correct dimensions
if size(time,1) == 1; time = time'; end
if size(depth,1) == 1; depth = depth'; end
if size(temp,1) == 1; temp = temp'; end
if size(lat,1) == 1; lat = lat'; end
if size(lon,1) == 1; lon = lon'; end
if size(ax,1) == 1; ax = ax'; end
if size(ay,1) == 1; ay = ay'; end
if size(az,1) == 1; az = az'; end
%if csv or xlsx files were used, break the larger matrix into variables
%based on format outlined in 'Introduction'
if exist('Data','var') %if the Data variable exist (xlsx or csv read in)
    %split it into variables
    time = Data(:,1); depth = Data(:,2); temp = Data(:,5); %time, depth, temperature
    lat = Data(:,3); lon = Data(:,4); %latitude and longitude
    ax = Data(:,6); ay = Data(:,7); az = Data(:,8); %acceleration
end

%Date and time information
%find day(s), calendar month(s), and calendar year(s) of each
%variable(depth, temp, lat/lon, etc.) measurement
timeYear = time.Year; %deployment year(s); same indices as 'time' vector for identifying year changes in variable data
timeRow = (mod(timeYear,4) ~= 0).*1 + (mod(timeYear,4) == 0).*2;  %determine row (1:Common 2:Leap) for all deployment measurements
timeMonth = time.Month; %deployment month(s); same indices as 'time' vector for identifying month changes in variable data
timeDay = monthStartDay(sub2ind(size(monthStartDay),timeRow,timeMonth)) + time.Day - 1; %deployment day(s) in day of year
%find day(s), calendar month(s), and calendar year(s) of deployment to pass to main model
Year = unique(timeYear,'stable'); %year(s) of deployment
Month = [timeMonth(diff(timeMonth) ~= 0); timeMonth(end)]; %month(s) of deployment
Day = timeDay(diff(timeDay) ~= 0); %days of deployment

%%%%%%if the data set resides within one calendar month
if length(Month) == 1
    %panel tilt and azimuth are not recorded directly and must be calculated
    pt = paneltilt(ax,ay,az); %panel tilt calculations
    pa = panelazimuth(lat, lon); %panel azimuth calculations
    %interpolate file contents
    Time = dateshift(time(1), 'start', 'day'):(1/86400):dateshift(time(end), 'end', 'day')-(1/86400); %full time vector, second resolution
    [Depth, Latitude, Longitude, T, Panel_Tilt, Panel_Azimuth] = cust_interp(time, Time, depth, lat, lon, temp, pt, pa);
    %write out a mat file with the necessary model inputs
    filename = sprintf('%s_model.mat',base); %pass to main model
    save(filename,'Time','Depth','Latitude','Longitude','T','Panel_Tilt','Panel_Azimuth','-mat'); %save model read matlab file
end

%%%%%%if data set spans multiple calendar months or years
%find the year start and end in the variable measurement indices to
%ensure strt and en indices only pull from current year
timeYear_strt = zeros(size(Year)); timeYear_end = zeros(size(Year)); %preallocate matrices
for iyear = 1:length(Year) %for each year deployment spans
    timeYear_strt(iyear) = find(timeYear == Year(iyear),1,'first');
    timeYear_end(iyear) = find(timeYear == Year(iyear),1,'last');
end
yearPointer = 1; %initialize 'Year' vector index pointer
filename = string(nan(length(Month),1)); %preallocating a matrix to save filenames
for iMonth = 1:length(Month) %iterate through the calendar months of the deployment data
    if iMonth >= 2 && Month(iMonth) == 1; yearPointer = yearPointer + 1; end %interate pointer if new year detected
    %extract month chunks using variable start and end indices
    if iMonth == 1 %if it is the first month
        strt = 1; %first measurement of deployment
        en = find(timeMonth(timeYear_strt(yearPointer):timeYear_end(yearPointer)) == Month(iMonth),1,'last'); %last measurement of first month
    elseif iMonth == length(Month) %if it is the last experienced month
        rel_timeMonth = timeMonth(timeYear_strt(yearPointer):timeYear_end(yearPointer)); %Month vector for only current year
        rel_time = time(timeYear_strt(yearPointer):timeYear_end(yearPointer)); %used to locate relative month value in full deployment
        rel_strt = find(rel_timeMonth == Month(iMonth),1,'first'); %first measurement of month
        strt = find(time == rel_time(rel_strt)); %first measurement of month in deployment indices
        en = length(time); %last measurement of deployment
    else %for full months
        rel_timeMonth = timeMonth(timeYear_strt(yearPointer):timeYear_end(yearPointer)); %Month vector for only current year
        rel_time = time(timeYear_strt(yearPointer):timeYear_end(yearPointer)); %used to locate relative month value in full deployment
        rel_strt = find(rel_timeMonth == Month(iMonth),1,'first'); %first measurement of month; relative to current year
        rel_en = find(rel_timeMonth == Month(iMonth),1,'last'); %last measurement of month; realtive to current year
        strt = find(time == rel_time(rel_strt)); %first measurement of month in deployment indices
        en = find(time ==rel_time(rel_en)); %last measurement of month in deployment indices
    end
    %seperate out the current month of data
    time_month = time(strt:en);
    depth_month = depth(strt:en);
    lat_month = lat(strt:en);
    lon_month = lon(strt:en);
    temp_month = temp(strt:en);
    ax_month = ax(strt:en);
    ay_month = ay(strt:en);
    az_month = az(strt:en);
    %calculate pt and pa
    pt_month = paneltilt(ax_month,ay_month,az_month); %panel tilt calculations
    pa_month = panelazimuth(lat_month, lon_month); %panel azimuth calculations
    %interpolate
    Time = dateshift(time_month(1), 'start', 'day'):(1/86400):dateshift(time_month(end), 'end', 'day')-(1/86400); %full time vector, second resolution
    [Depth, Latitude, Longitude, T, Panel_Tilt, Panel_Azimuth] = cust_interp(time_month, Time, depth_month, lat_month, lon_month, temp_month, pt_month, pa_month);
    %write to file
    %pad single numbered months with leading zero
    if length(num2str(Month(iMonth))) == 1; monthPrint = sprintf('0%d',Month(iMonth));
    else; monthPrint = string(Month(iMonth)); end
    filename(iMonth) = sprintf('%s_%d_%s_model.mat',base,Year(yearPointer),monthPrint); %filenames to pass to main model
    save(filename(iMonth),'Time','Depth','Latitude','Longitude','T','Panel_Tilt','Panel_Azimuth','-mat')
end

% Subfunction Definitions
%Panel Tilt Calculations
function pt = paneltilt(ax,ay,az)
%{
To calculate the panel tilt from the acceleration data, we need to find
the points where the proper acceleration (measured accel/9.81) is equal to
one; suggesting that the acceleration experienced by the device was only
caused by earth's gravity (9.81/9.81 = 1). However, if the strict condition
of proper acceleration is equal to one is applied to the data, very few
data points are left to interpolate. Therefore, we introduced a tolerance
to the proper acceleration limit (1 + aTol) and do a cost analysis of
the amount of data a tolerance adds vs the angle error it causes. The
maximum angle error is used in this comparison, which occurs when the
acceleration tolerance is added to the proper acceleration orthogonally.
%}
prop_a = sqrt(ax.^2 + ay.^2 + az.^2)./9.81; %magnitude of proper acceleration from data
aTol = 0:0.01:1; %range of tolerances for proper acceleration limit
%Cost analysis of tolerance error vs data inclusion (percent exclusion is
%used here for plotting/visualiztion purposes)
p_exc = zeros(size(aTol));    %preallocation of percent excluded vector
for i = 1:length(aTol)    %find percent of excluded data for each tolerance
    n_exc = sum(prop_a > (1+aTol(i)) | prop_a < (1-aTol(i))); %number of excluded points
    p_exc(i) = n_exc*100/length(prop_a); %percent of excluded points
end
inc95 = find(p_exc < 5,1,'first');  %find 95% inclusion marker
aTol95 = aTol(inc95);

%{
Assuming the default value of 95% inclusion of the data points does not
produce an unacceptable error, the following process will calculate the
panel tilt based on the valid acceleration data. (Valid is refering to data
points whose proper acceleration is within tolerance as described above)
%}
n = [0 0 -1]; %gravity vector of undisturbed accelerometer in Universal coordinate
ivalid_a = find(prop_a < (1+aTol95) & prop_a > (1-aTol95)); %find indices of valid acceleration measurements (1+/- tolerance)
pt = nan(size(prop_a)); %preallocation of panel tilt (degrees)
h = waitbar(0,'Calculating panel tilt...'); %create progress bar
for i = 1:length(ivalid_a) %calculate panel tilt at every valid acceleration measurement
    a = [-ay(ivalid_a(i)), ax(ivalid_a(i)), az(ivalid_a(i))]; %translate measured acceleration coordinates to aircraft coordinate system
    pt(ivalid_a(i)) = acosd(dot(n,a)/(norm(n)*norm(a))); %calculate panel tilt using dot product definition
    pcomp = i/length(ivalid_a); %percent complete
    waitbar(pcomp) %increment progress bar
end
close(h) %close progress bar
end

%Panel Azimuth Calculations
function pa = panelazimuth(lat, lon)
%{
To calculate the panel azimuth, we assume linear movment between latitude
and longitude measurements. Assume lon 0 to 360
%}
%find lat lon direction travelled, used for classifying quadrant
lat_diff = diff(lat);
lon_diff = diff(lon);
%find km equivalent of lat lon diff
lat_km = diff(deg2km(lat,'earth'));
lon_km = diff(deg2km(lon,'earth'));
h = waitbar(0,'Calculating panel azimuth...'); %create progress bar

adj = abs(lat_km); %distance adjacent to theta angle
opp = abs(lon_km); %distance opposite to theta angle
theta = atand(opp./adj);  %angle of panel relative to quadrant
pa = 0.*(lat_diff > 0 & lon_diff == 0) + ... %moves directly N
    180.*(lat_diff < 0 & lon_diff == 0) + ... %moves directly S
    90.*(lat_diff == 0 & lon_diff > 0) + ... %moves directly W
    270.*(lat_diff == 0 & lon_diff < 0) + ... %moves directly E
    (360-theta).*(lat_diff > 0 & lon_diff > 0) + ... %second quadrant
    theta.*(lat_diff > 0 & lon_diff < 0) + ... %first quadrant
    (180+theta).*(lat_diff < 0 & lon_diff > 0) + ... %third quadrant
    (180-theta).*(lat_diff < 0 & lon_diff < 0); %fourth quadrant
pa = [pa(1); pa];
close(h) %close progress bar
end

%Interpolation Function
function [Depth, Latitude, Longitude, T, Panel_Tilt, Panel_Azimuth] = cust_interp(time, Time, depth, lat, lon, temp, pt, pa)
%ensure unique time points
[time,ui,ti] = unique(time);
%create vectors to correspond with unique time vector
depth = depth(ui);
lat = lat(ui);
lon = lon(ui);
temp = temp(ui);
pt = pt(ui);
pa = pa(ui);
%interpolate out to model resolution
Depth = interp1(time, depth, Time);
Latitude = interp1(time, lat, Time);
Longitude = interp1(time, lon, Time);
T = interp1(time, temp, Time);
Panel_Tilt = interp1(time, pt, Time);
Panel_Azimuth = interp1(time, pa, Time);
end

end