function [Zenith, Solar_Azimuth, Incidence_Angle] = ...
    Solar_Position(latitude,longitude,datetime,Panel_Tilt_Bound, Panel_Azimuth,varargin)
%% Solar Position Algorithm (SPA)
%PROGRAMMER: Collin Krawczyk
%DATE: 5/20/2019
%DETAILS: Calculates zenith and azimuth angles given a location
%and datetime.
%Algorithm taken from https://www.nrel.gov/docs/fy08osti/34302.pdf
%Required Inputs: 
%   latitude            The latitude of interest with any amount of decimal
%                       places. Positive or negative if north or south of
%                       the equator, respectively. Limits from -90 to 90
%                       degrees.
%   longitude           The longitude of interest with any amount of decimal
%                       places. Positive or negative for east or west of 
%                       Greenwich, respectively. Limits from -180 to 180
%                       degrees.
%   datetime            Datetime format that includes year, month, day,
%                       hour, minute, and second. 
%   Panel_Tilt_Bound    The bound on the panel tilt in degrees.
%   Panel_Azimuth       The panel azimuth in degrees.
%Optional Input:
%   elevation           The elevation of interest in meters. Can be left 
%                       blank for sea level calculation (0 meters) and at
%                       standard temperature and pressure (0 Celcius, 1
%                       atm). If elevation is inputted, must need aerospace
%                       toolbox to get temperature and pressure at the
%                       elevation.
%Outputs:
%   Azimuth             The solar azimuth in degrees. Measured eastward
%                       from north. Accurate to 0.0001 degrees if compared 
%                       to SPA example. If compared to NOAA solar 
%                       calculation, it is accurate to 0.001 degrees.
%                       Issue comes from calculations of deltaT 
%                       (equation 2 on page 2). If a more precise method of 
%                       calculating deltaT can be implemented, the error
%                       can be reduced further. The example in SPA had a 
%                       deltaT of 67s. The calcuated deltaT was 64.5078s 
%                       using online NASA equations. 
%   Zenith              The solar zenith in degrees. Measured from the
%                       vertical to the sun's location. Accurate to 0.03 
%                       degrees if compared to SPA example due to issues 
%                       with varying pressure, temperature, and elevation. 
%                       If compared to NOAA solar calculation, it is 
%                       accurate to 0.05 degrees.
%   Incidence_Angle     The angle between a ray incident on the panel and
%                       the perpendicular surface in degrees

%% Input Variables

%Parse the elevation variable

p = inputParser;

default_ele = 0; %0 meters of elevation (sea level)

%Required Inputs
addRequired(p,'latitude',@isnumeric);
addRequired(p,'longitude',@isnumeric);
addRequired(p,'datetime', @isdatetime)

%Optional Input
addOptional(p,'elevation',default_ele,@isnumeric);

parse(p,latitude,longitude,datetime,varargin{:});

%Parser degugging check
% if ~isempty(fieldnames(p.Unmatched))
%     disp('Extra inputs: ');
%     disp(p.Unmatched)
% end
% if ~isempty(p.UsingDefaults)
%     disp('Using Defaults: ')
%     disp(p.UsingDefaults)
% end

Inputs.elevation = p.Results.elevation;

%Set the other variables equal to variables inputted
Inputs.datetime = datetime;
[Inputs.year,Inputs.month,Inputs.day] = ymd(Inputs.datetime);
[Inputs.hour,Inputs.minute,Inputs.second] = hms(Inputs.datetime);
Inputs.latitude = latitude; %Latitude from -90 to 90 degrees
Inputs.longitude = longitude; %Longitude from -180 to 180 degrees

%% Convert each input variable into equal length arrays

%Year

if length(Inputs.year) > length(Inputs.month) || ...
        length(Inputs.year) > length(Inputs.day) || ...
        length(Inputs.year) > length(Inputs.hour) || ...
        length(Inputs.year) > length(Inputs.minute) || ...
        length(Inputs.year) > length(Inputs.second) || ...
        length(Inputs.year) > length(Inputs.latitude) || ...
        length(Inputs.year) > length(Inputs.longitude) || ...
        length(Inputs.year) > length(Inputs.elevation)
    
    Inputs.month = repmat(Inputs.month(1),length(Inputs.year),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.year),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.year),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.year),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.year),1)';
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.year),1)';
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.year),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.year),1)';
    
end

%Month

if length(Inputs.month) > length(Inputs.year) || ...
        length(Inputs.month) > length(Inputs.day) || ...
        length(Inputs.month) > length(Inputs.hour) || ...
        length(Inputs.month) > length(Inputs.minute) || ...
        length(Inputs.month) > length(Inputs.second) || ...
        length(Inputs.month) > length(Inputs.latitude) || ...
        length(Inputs.month) > length(Inputs.longitude) || ...
        length(Inputs.month) > length(Inputs.elevation)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.month),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.month),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.month),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.month),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.month),1)';
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.month),1)';
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.month),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.month),1)';
end

%Day

if length(Inputs.day) > length(Inputs.year) || ...
        length(Inputs.day) > length(Inputs.month) || ...
        length(Inputs.day) > length(Inputs.hour) || ...
        length(Inputs.day) > length(Inputs.minute) || ...
        length(Inputs.day) > length(Inputs.second) || ...
        length(Inputs.day) > length(Inputs.latitude) || ...
        length(Inputs.day) > length(Inputs.longitude) || ...
        length(Inputs.day) > length(Inputs.elevation)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.day),1)';
    Inputs.month = repmat(Inputs.month(1),length(Inputs.day),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.day),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.day),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.day),1)';
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.day),1)';
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.day),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.day),1)';
end

%Hour

if length(Inputs.hour) > length(Inputs.year) || ...
        length(Inputs.hour) > length(Inputs.month) || ...
        length(Inputs.hour) > length(Inputs.day) || ...
        length(Inputs.hour) > length(Inputs.minute) || ...
        length(Inputs.hour) > length(Inputs.second) || ...
        length(Inputs.hour) > length(Inputs.latitude) || ...
        length(Inputs.hour) > length(Inputs.longitude) || ...
        length(Inputs.hour) > length(Inputs.elevation)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.hour),1)';
    Inputs.month = repmat(Inputs.month(1),length(Inputs.hour),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.hour),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.hour),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.hour),1)';
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.hour),1)';
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.hour),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.hour),1)';
end

%Minute

if length(Inputs.minute) > length(Inputs.year) || ...
        length(Inputs.minute) > length(Inputs.month) || ...
        length(Inputs.minute) > length(Inputs.day) || ...
        length(Inputs.minute) > length(Inputs.hour) || ...
        length(Inputs.minute) > length(Inputs.second) || ...
        length(Inputs.minute) > length(Inputs.latitude) || ...
        length(Inputs.minute) > length(Inputs.longitude) || ...
        length(Inputs.minute) > length(Inputs.elevation)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.minute),1)';
    Inputs.month = repmat(Inputs.month(1),length(Inputs.minute),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.minute),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.minute),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.minute),1)';
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.minute),1)';
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.minute),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.minute),1)';
end

%Second

if length(Inputs.second) > length(Inputs.year) || ...
        length(Inputs.second) > length(Inputs.month) || ...
        length(Inputs.second) > length(Inputs.day) || ...
        length(Inputs.second) > length(Inputs.hour) || ...
        length(Inputs.second) > length(Inputs.minute) || ...
        length(Inputs.second) > length(Inputs.latitude) || ...
        length(Inputs.second) > length(Inputs.longitude) || ...
        length(Inputs.second) > length(Inputs.elevation)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.second),1)';
    Inputs.month = repmat(Inputs.month(1),length(Inputs.second),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.second),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.second),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.second),1)';
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.second),1)';
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.second),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.second),1)';
end

%Latitude

if length(Inputs.latitude) > length(Inputs.year) || ...
        length(Inputs.latitude) > length(Inputs.month) || ...
        length(Inputs.latitude) > length(Inputs.day) || ...
        length(Inputs.latitude) > length(Inputs.hour) || ...
        length(Inputs.latitude) > length(Inputs.minute) || ...
        length(Inputs.latitude) > length(Inputs.second) || ...
        length(Inputs.latitude) > length(Inputs.elevation)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.latitude),1)';
    Inputs.month = repmat(Inputs.month(1),length(Inputs.latitude),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.latitude),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.latitude),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.latitude),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.latitude),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.latitude),1)';
end
    
if length(Inputs.latitude) > length(Inputs.longitude)
    Inputs.latitude = Inputs.latitude;
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.latitude),1)';
end
    
if length(Inputs.latitude) == length(Inputs.longitude)
    Inputs.latitude = Inputs.latitude;
    Inputs.longitude = Inputs.longitude;
end
    
if length(Inputs.longitude) > length(Inputs.latitude)
    Inputs.longitude = Inputs.longitude;
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.longitude),1)';   
end


%Longitude

if length(Inputs.longitude) > length(Inputs.year) || ...
        length(Inputs.longitude) > length(Inputs.month) || ...
        length(Inputs.longitude) > length(Inputs.day) || ...
        length(Inputs.longitude) > length(Inputs.hour) || ...
        length(Inputs.longitude) > length(Inputs.minute) || ...
        length(Inputs.longitude) > length(Inputs.second) || ...
        length(Inputs.longitude) > length(Inputs.elevation)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.longitude),1)';
    Inputs.month = repmat(Inputs.month(1),length(Inputs.longitude),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.longitude),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.longitude),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.longitude),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.longitude),1)';
    Inputs.elevation = repmat(Inputs.elevation(1),length(Inputs.longitude),1)';
end

%Elevation

if length(Inputs.elevation) > length(Inputs.year) || ...
        length(Inputs.elevation) > length(Inputs.month) || ...
        length(Inputs.elevation) > length(Inputs.day) || ...
        length(Inputs.elevation) > length(Inputs.hour) || ...
        length(Inputs.elevation) > length(Inputs.minute) || ...
        length(Inputs.elevation) > length(Inputs.second) || ...
        length(Inputs.elevation) > length(Inputs.latitude) || ...
        length(Inputs.elevation) > length(Inputs.longitude)
    
    Inputs.year = repmat(Inputs.year(1),length(Inputs.elevation),1)';
    Inputs.month = repmat(Inputs.month(1),length(Inputs.elevation),1)';
    Inputs.day = repmat(Inputs.day(1),length(Inputs.elevation),1)';
    Inputs.hour = repmat(Inputs.hour(1),length(Inputs.elevation),1)';
    Inputs.minute = repmat(Inputs.minute(1),length(Inputs.elevation),1)';
    Inputs.second = repmat(Inputs.second(1),length(Inputs.elevation),1)';
    Inputs.latitude = repmat(Inputs.latitude(1),length(Inputs.elevation),1)';
    Inputs.longitude = repmat(Inputs.longitude(1),length(Inputs.elevation),1)';
end

%% Check to make sure inputs are between ranges

for i = 1:length(Inputs.year)
    if Inputs.year(i) < -2000 || Inputs.year(i) > 6000
        warning('Results may be inaccurate. The SPA is valid for years -2000 to 6000');
    end

    if Inputs.month(i) < 1 || Inputs.month(i) > 12
        error('Error:outsideBoundaries',...
        '\nThe month must be between 1 and 12 where 1 corresponds to January, 2 corresponds to February, etc.');
    end

    if Inputs.day(i) < 1 || Inputs.day(i) > 31
        error('Error:outsideBoundaries',...
        '\nThe day must be between 1 and 31 where 1 corresponds to day 1, 2 corresponds to day 2, etc.');
    end

    if Inputs.hour(i) < 0 || Inputs.hour(i) > 24
        error('Error:outsideBoundaries',...
        '\nThe hour must be between 0 and 24 (UTC) where 0 corresponds to hour 0, 1 corresponds to hour 1, etc.');
    end

    if Inputs.minute(i) < 0 || Inputs.minute(i) > 60
        error('Error:outsideBoundaries',...
        '\nThe minute must be between 0 and 60 where 0 corresponds to 0 minutes, 1 corresponds to 1 minute, etc.');
    end

    if Inputs.second(i) < 0 || Inputs.second(i) > 60
        error('Error:outsideBoundaries',...
        '\nThe second must be between 0 and 60 where 0 corresponds to 0 seconds, 1 corresponds to 1 second, etc.');
    end

    if Inputs.latitude(i) < -90 || Inputs.latitude(i) > 90
        error('Error:outsideBoundaries',...
        '\nThe latitude must be between -90 and 90 degrees.');
    end

    if Inputs.longitude(i) < -180 || Inputs.longitude(i) > 180
        error('Error:outsideBoundaries',...
        '\nThe longitude must be between -180 and 180 degrees.');
    end
end
  
%% Calculate pressure and temperature based on elevation
%Need aerospace toolbox for atmosisa function

Temperature = zeros(length(Inputs.elevation),1);
Pressure = zeros(length(Inputs.elevation),1);

for i = 1:length(Inputs.elevation)
    
    if Inputs.elevation(i) == 0
        Inputs.pressure(i) = 101325/100; %Standard Pressure in millibars
        Inputs.temperature(i) = 0; %Standard Temperature in Celcius
        
    else
        [Temperature(i),Inputs.sound_speed(i),Pressure(i),Inputs.density(i)]...
            = atmosisa(Inputs.elevation(i)); %Aerospace Toolbox
    
        Inputs.pressure(i) = Pressure(i)/100; %Convert to millibar
        Inputs.temperature(i) = Temperature(i)-273.15; %Convert to Celcius
    end
end

%% CALCULATE THE JULIAN DAY (JD) using NREL equation (4)

%Day of the month with decimal time 
%(e.g. for the second day of the month at 12:30:30 UT, D = 2.521180556);

date.decimal_day = Inputs.day+(Inputs.hour/24)+ ...
    (Inputs.minute*(0.04166666666666666/60))+ ...
    (Inputs.second*(0.000694444444444444/60)); 

for i = 1:length(Inputs.month)
    if (Inputs.month(i) < 3)
        Inputs.month(i) = Inputs.month(i) + 12;
        Inputs.year(i) = Inputs.year(i) - 1;
    end
end

date.jd = fix(365.25.*(Inputs.year+4716)) + fix(30.6001.* ...
    (Inputs.month+1))+date.decimal_day-1524.5;

p = zeros(1,length(date.jd));
for i = 1:length(date.jd)
    if date.jd(i) > 2299160
        p(i) = fix(Inputs.year(i)/100);
        date.jd(i) = date.jd(i) + (2-p(i)+fix(p(i)/4));
    end
end

%% CALCULATE DIFFERENCE BETWEEN UTC AND TERRESTRIAL TIME (DELTA T IN NREL)
%(will be in seconds) 
%All taken from https://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html

%Year decimal that gives year decimal for the middle of the month

for i = 1:length(Inputs.year)
    date.decimal_year(i,1) = Inputs.year(i) + (Inputs.month(i)-0.5)/12;
end

j = zeros(length(Inputs.year),1);
for i = 1:length(Inputs.year)
    if Inputs.year(i) < -500 || Inputs.year(i) > 2150
        j(i) = (Inputs.year(i)-1820)/100;
        date.time_difference(i,1) = -20 + 32*j(i).^2;
        
    elseif Inputs.year(i) >= -500 && Inputs.year(i) < 500
        j(i) = date.decimal_year(i)/100;
        date.time_difference(i,1) = 10583.6 - 1014.41*j(i) + ...
            33.78311*j(i).^2 - 5.952053*j(i).^3 - 0.1798452*j(i).^4 + ...
            0.022174192*j(i).^5 + 0.0090316521*j(i).^6;
        
    elseif Inputs.year(i) >= 500 && Inputs.year(i) < 1600
        j(i) = (date.decimal_year(i)-1000)/100;
        date.time_difference(i,1) = 1574.2 - 556.01*j(i) + ...
            71.23472*j(i).^2 + 0.319781*j(i).^3 - 0.8503463*j(i).^4 - ...
            0.005050998*j(i).^5 + 0.0083572073*j(i).^6;
        
    elseif Inputs.year(i) >= 1600 && Inputs.year(i) < 1700
        j(i) = date.decimal_year(i) - 1600;
        date.time_difference(i,1) = 120 - 0.9808*j(i) - 0.01532*j(i).^2 ...
            + (j(i).^3)/7129;
        
    elseif Inputs.year(i) >= 1700 && Inputs.year(i) < 1800
        j(i) = date.decimal_year(i) - 1700;
        date.time_difference(i,1) = 8.83 + 0.1603*j(i) - ...
            0.0059285*j(i).^2 + 0.00013336*j(i).^3 - (j(i).^4)/1174000;
        
    elseif Inputs.year(i) >= 1800 && Inputs.year(i) < 1860
        j(i) = date.decimal_year(i) - 1800;
        date.time_difference(i,1) = 13.72 - 0.332447*j(i) + ...
            0.0068612*j(i).^2 + 0.0041116*j(i).^3 - 0.00037436*j(i).^4 ...
            + 0.0000121272*j(i).^5 - 0.0000001699*j(i).^6 + ...
            0.000000000875*j(i).^7;
        
    elseif Inputs.year(i) >= 1860 && Inputs.year(i) < 1900
        j(i) = date.decimal_year(i) - 1860;
        date.time_difference(i,1) = 7.62 + 0.5737*j(i) - ...
            0.251754*j(i).^2 + 0.01680668*j(i).^3 - ...
            0.0004473624*j(i).^4 + (j(i).^5)/233174;
        
    elseif Inputs.year(i) >= 1900 && Inputs.year(i) < 1920
        j(i) = date.decimal_year(i) - 1900;
        date.time_difference(i,1) = -2.79 + 1.494119*j(i) - ...
            0.0598939*j(i).^2 + 0.0061966*j(i).^3 - 0.000197*j(i).^4;
        
    elseif Inputs.year(i) >= 1920 && Inputs.year(i) < 1941
        j(i) = date.decimal_year(i) - 1920;
        date.time_difference(i,1) = 21.20 + 0.84493*j(i) - ...
            0.076100*j(i).^2 + 0.0020936*j(i).^3;
        
    elseif Inputs.year(i) >= 1941 && Inputs.year(i) < 1961
        j(i) = date.decimal_year(i) - 1950;
        date.time_difference(i,1) = 29.07 + 0.407*j(i) - ...
            (j(i).^2)/233 + (j(i).^3)/2547;
        
    elseif Inputs.year(i) >= 1961 && Inputs.year(i) < 1986
        j(i) = date.decimal_year(i) - 1975;
        date.time_difference(i,1) = 45.45 + 1.067*j(i) - ...
            (j(i).^2)/260 - (j(i).^3)/718;
        
    elseif Inputs.year(i) >= 1986 && Inputs.year(i) < 2005
        j(i) = date.decimal_year(i) - 2000;
        date.time_difference(i,1) = 63.86 + 0.3345*j(i) - ...
            0.060374*j(i).^2 + 0.0017275*j(i).^3 + 0.000651814*j(i).^4 ...
            + 0.00002373599*j(i).^5;
        
    elseif Inputs.year(i) >= 2005 && Inputs.year(i) < 2050
        j(i) = date.decimal_year(i) - 2000;
        date.time_difference(i,1) = 62.92 + 0.32217*j(i) + ...
            0.005589*j(i).^2;
        
    elseif Inputs.year(i) >= 2050 && Inputs.year(i) <= 2150
        date.time_difference(i,1) = -20 + ...
            32*((date.decimal_year(i)-1820)/100).^2 - ...
            0.5628*(2150-date.decimal_year(i));
    end
end

%% CALCULATE NECESSARY JULIAN DAY CALCULATIONS 

%Julian Ephemeris Day (JDE) using NREL equation (5)

date.JDE = date.jd+(date.time_difference'/86400); 

%Julian Century (JC) using NREL equation (6)

date.JC = (date.jd-2451545)/36525; 

%Julian Ephemeris Century (JCE) using NREL equation (7)

date.JCE = (date.JDE-2451545)/36525; 

%Julian Ephemeris Millennium (JME) for 2000 standard epoch 
%using NREL equation (8)

date.JME = (date.JCE)/10; 

%% INPUT CONSTANTS A, B, AND C FROM TABLE A4.2 IN NREL

%Load in the .mat file that corresponds to Table A4.2
%This could be localized by inputting these into this script instead of
%having a seperate .mat file

load('TableA4.2.mat');

%                               L ARRAYS
%Create L0 arrays

L.L0.L0A = zeros(64,1);
L.L0.L0B = zeros(64,1);
L.L0.L0C = zeros(64,1);
for i = 1:64
    L.L0.L0A(i) = TableA4_2{1,2}((i),1);
    L.L0.L0B(i) = TableA4_2{1,3}((i),1);
    L.L0.L0C(i) = TableA4_2{1,4}((i),1);
end

%CreatE L1 arrays

L.L1.L1A = zeros(34,1);
L.L1.L1B = zeros(34,1);
L.L1.L1C = zeros(34,1);
for i = 65:98
    L.L1.L1A(i-64) = TableA4_2{1,2}((i),1);
    L.L1.L1B(i-64) = TableA4_2{1,3}((i),1);
    L.L1.L1C(i-64) = TableA4_2{1,4}((i),1);
end

%Create L2 arrays

L.L2.L2A = zeros(20,1);
L.L2.L2B = zeros(20,1);
L.L2.L2C = zeros(20,1);
for i = 99:118
    L.L2.L2A(i-98) = TableA4_2{1,2}((i),1);
    L.L2.L2B(i-98) = TableA4_2{1,3}((i),1);
    L.L2.L2C(i-98) = TableA4_2{1,4}((i),1);
end

%Create L3 arrays

L.L3.L3A = zeros(7,1);
L.L3.L3B = zeros(7,1);
L.L3.L3C = zeros(7,1);
for i = 119:125
    L.L3.L3A(i-118) = TableA4_2{1,2}((i),1);
    L.L3.L3B(i-118) = TableA4_2{1,3}((i),1);
    L.L3.L3C(i-118) = TableA4_2{1,4}((i),1);
end

%Create L4 arrays

L.L4.L4A = zeros(3,1);
L.L4.L4B = zeros(3,1);
L.L4.L4C = zeros(3,1);
for i = 126:128
    L.L4.L4A(i-125) = TableA4_2{1,2}((i),1);
    L.L4.L4B(i-125) = TableA4_2{1,3}((i),1);
    L.L4.L4C(i-125) = TableA4_2{1,4}((i),1);
end

%Create L5 arrays

L.L5.L5A = TableA4_2{1,2}((129),1);
L.L5.L5B = TableA4_2{1,3}((129),1);
L.L5.L5C = TableA4_2{1,4}((129),1);

%                             B ARRAYS
%Create B0 arrays

B.B0.B0A = zeros(5,1);
B.B0.B0B = zeros(5,1);
B.B0.B0C = zeros(5,1);
for i = 130:134
    B.B0.B0A(i-129) = TableA4_2{1,2}((i),1);
    B.B0.B0B(i-129) = TableA4_2{1,3}((i),1);
    B.B0.B0C(i-129) = TableA4_2{1,4}((i),1);
end

%Create B1 arrays

B.B1.B1A = zeros(2,1);
B.B1.B1B = zeros(2,1);
B.B1.B1C = zeros(2,1);
for i = 135:136
    B.B1.B1A(i-134) = TableA4_2{1,2}((i),1);
    B.B1.B1B(i-134) = TableA4_2{1,3}((i),1);
    B.B1.B1C(i-134) = TableA4_2{1,4}((i),1);
end

%                             R ARRAYS
%Create R0 arrays

R.R0.R0A = zeros(40,1);
R.R0.R0B = zeros(40,1);
R.R0.R0C = zeros(40,1);
for i = 137:176
    R.R0.R0A(i-136) = TableA4_2{1,2}((i),1);
    R.R0.R0B(i-136) = TableA4_2{1,3}((i),1);
    R.R0.R0C(i-136) = TableA4_2{1,4}((i),1);
end

%Create R1 arrays

R.R1.R1A = zeros(10,1);
R.R1.R1B = zeros(10,1);
R.R1.R1C = zeros(10,1);
for i = 177:186
    R.R1.R1A(i-176) = TableA4_2{1,2}((i),1);
    R.R1.R1B(i-176) = TableA4_2{1,3}((i),1);
    R.R1.R1C(i-176) = TableA4_2{1,4}((i),1);
end

%Create R2 arrays

R.R2.R2A = zeros(6,1);
R.R2.R2B = zeros(6,1);
R.R2.R2C = zeros(6,1);
for i = 187:192
    R.R2.R2A(i-186) = TableA4_2{1,2}((i),1);
    R.R2.R2B(i-186) = TableA4_2{1,3}((i),1);
    R.R2.R2C(i-186) = TableA4_2{1,4}((i),1);
end

%Create R3 arrays

R.R3.R3A = zeros(2,1);
R.R3.R3B = zeros(2,1);
R.R3.R3C = zeros(2,1);
for i = 193:194
    R.R3.R3A(i-192) = TableA4_2{1,2}((i),1);
    R.R3.R3B(i-192) = TableA4_2{1,3}((i),1);
    R.R3.R3C(i-192) = TableA4_2{1,4}((i),1);
end

%Create R4 arrays

R.R4.R4A = TableA4_2{1,2}((195),1);
R.R4.R4B = TableA4_2{1,3}((195),1);
R.R4.R4C = TableA4_2{1,4}((195),1);

%% CALCULATE L'S, B'S, AND R'S USING NREL EQUATIONS (9) AND (10)

%                                 L's
%L0 

L.L0.ind = L.L0.L0A .* cos(L.L0.L0B + L.L0.L0C .* ...
            date.JME);
        
L.L0.sum = sum(L.L0.ind,1);

%L1

L.L1.ind = L.L1.L1A .* cos(L.L1.L1B + L.L1.L1C .* ...
            date.JME);
        
L.L1.sum = sum(L.L1.ind,1);

%L2

L.L2.ind = L.L2.L2A .* cos(L.L2.L2B + L.L2.L2C .* ...
            date.JME);
        
L.L2.sum = sum(L.L2.ind,1);

%L3

L.L3.ind = L.L3.L3A .* cos(L.L3.L3B + L.L3.L3C .* ...
            date.JME);
        
L.L3.sum = sum(L.L3.ind,1);

%L4

L.L4.ind = L.L4.L4A .* cos(L.L4.L4B + L.L4.L4C .* ...
            date.JME);
        
L.L4.sum = sum(L.L4.ind,1);
        
%L5

L.L5.ind = L.L5.L5A .* cos(L.L5.L5B + L.L5.L5C .* ...
            date.JME);
        
L.L5.sum = sum(L.L5.ind,1);

%                                   B's
%B0

B.B0.ind = B.B0.B0A .* cos(B.B0.B0B + B.B0.B0C .* ...
            date.JME);
        
B.B0.sum = sum(B.B0.ind,1);

%B1

B.B1.ind = B.B1.B1A .* cos(B.B1.B1B + B.B1.B1C .* ...
            date.JME);
        
B.B1.sum = sum(B.B1.ind,1);

%                               R's
%R0

R.R0.ind = R.R0.R0A .* cos(R.R0.R0B + R.R0.R0C .* ...
            date.JME);
        
R.R0.sum = sum(R.R0.ind,1);

%R1

R.R1.ind = R.R1.R1A .* cos(R.R1.R1B + R.R1.R1C .* ...
            date.JME);
        
R.R1.sum = sum(R.R1.ind,1);

%R2

R.R2.ind = R.R2.R2A .* cos(R.R2.R2B + R.R2.R2C .* ...
            date.JME);
        
R.R2.sum = sum(R.R2.ind,1);

%R3

R.R3.ind = R.R3.R3A .* cos(R.R3.R3B + R.R3.R3C .* ...
            date.JME);
        
R.R3.sum = sum(R.R3.ind,1);

%R4

R.R4.ind = R.R4.R4A .* cos(R.R4.R4B + R.R4.R4C .* ...
            date.JME);
        
R.R4.sum = sum(R.R4.ind,1);

%% CALCULATE TOTAL L, B, AND R USING NREL EQUATIONS (11) AND (12)

%L in radians using NREL section 3.2.4

L.radians =   (L.L0.sum + L.L1.sum .* (date.JME) + ...
    L.L2.sum .* (date.JME).^2 + L.L3.sum .* (date.JME).^3 + ...
    L.L4.sum .* (date.JME).^4 + L.L5.sum .* (date.JME).^5)./10^8;

%L in degrees using NREL section 3.2.5

L.degrees = (L.radians*180)/pi;

%L normalized to 360 degrees using NREL section 3.2.6

limit.F = abs(rem(L.degrees/360,1));

for i = 1:length(L.degrees)
    if L.degrees(i) >= 0
        L.normalized(i) = 360*limit.F(i);
    elseif L.degrees(i) < 0 
        L.normalized(i) = 360-360*limit.F(i);
    end
end
    
%B using NREL section 3.2.7

B.radians = (B.B0.sum + B.B1.sum .*date.JME)./10^8;

B.degrees = (B.radians*180)/pi; 

%R using NREL section 3.2.8 (AU)

R.AU = (R.R0.sum + R.R1.sum .* (date.JME) + ...
    R.R2.sum .* (date.JME).^2 + ...
    R.R3.sum .* (date.JME).^3 + ...
    R.R4.sum .* (date.JME).^4)./10^8;

%% CALCULATE GEOCENTRIC LONGITUDE AND LATITUDE (THETA AND BETA)

%Theta (degrees) using NREL equation (13)

geocentric.longitude = L.normalized + 180;

%Limit to 360 degrees

limit.K = abs(rem(geocentric.longitude/360,1));

for i = 1:length(geocentric.longitude)
    if geocentric.longitude(i) >= 0
        geocentric.normalized.longitude(i) = 360*limit.K(i);
    elseif geocentric.longitude_nrel(i) < 0 
        geocentric.normalized.longitude(i) = 360-360*limit.K(i);
    end
end

%Beta (degrees);

geocentric.latitude = -1*B.degrees;

%% CALCULATE NUTATION IN LONGITUDE AND OBLIQUITY (DELTA PSI & DELTA EPSILON)

%Mean elongation of the moon from the sun (degrees) NREL equation (15)

X.X0 = 297.85036 + 445267.111480 .* (date.JCE) - 0.0019142 .* ...
    (date.JCE).^2 + ((date.JCE).^3) ./ 189474;

%Mean anomaly of the sun (Earth)(degrees) NREL equation (16)

X.X1 = 357.52772 + 35999.050340 .* (date.JCE) - 0.0001603 .* ...
    (date.JCE).^2 - ((date.JCE).^3) ./ 300000;

%Mean anomoly of the moon (degrees) NREL equation (17)

X.X2 = 134.96298 + 477198.867398 .* (date.JCE) + 0.0086972 .* ...
    (date.JCE).^2 + ((date.JCE).^3) ./ 56250;

%Moon argument of latitude (degrees) NREL equation (18)

X.X3 = 93.27191 + 483202.017538 .* (date.JCE) - 0.0036825 .* ...
    (date.JCE).^2 + ((date.JCE).^3) ./ 327270;

%Longitude of the assending node of the moon's mean orbit on the ecliptic,
%measured from mean equinox of the date (degrees) NREL equation (19)

X.X4 = 125.04452 - 1934.136261 .* (date.JCE) + 0.0020708 .* ...
    (date.JCE).^2 + ((date.JCE).^3) ./ 450000;

%% Load in Table A4.3 and calculate delta PSI and delta EPSILON

%Input Y's into a .mat file
% Y0 = xlsread('Solar Position.xlsx', 'Table 5', 'A3:A65');
% Y1 = xlsread('Solar Position.xlsx', 'Table 5', 'B3:B65');
% Y2 = xlsread('Solar Position.xlsx', 'Table 5', 'C3:C65');
% Y3 = xlsread('Solar Position.xlsx', 'Table 5', 'D3:D65');
% Y4 = xlsread('Solar Position.xlsx', 'Table 5', 'E3:E65');
% a = xlsread('Solar Position.xlsx', 'Table 5', 'F3:F65');
% b = xlsread('Solar Position.xlsx', 'Table 5', 'G3:G65');
% c = xlsread('Solar Position.xlsx', 'Table 5', 'H3:H65');
% d = xlsread('Solar Position.xlsx', 'Table 5', 'I3:I65');
% TableA4_3 = [Y0,Y1,Y2,Y3,Y4,a,b,c,d];
% save('TableA4.3.mat', 'TableA4_3');

load('TableA4.3.mat');

%Calculate delta PSI and delta EPSILON terms of 
%NREL equations (20) and (21)

C = struct2cell(X);

trig.terms = zeros(63,5,length(date.JME));

% Calculate the xy terms first

for i = 1:63
    for k = 1:5
        for j = 1:length(date.JCE)
            trig.terms(i,k,j)= C{k,1}(j) .* TableA4_3(i,k);
        end
    end
end

%Sum up xy terms and calculate psi and epsilon

for i = 1:63
     for j = 1:length(date.JCE)
         psi.delta.sum(i,1,j) = sum(trig.terms(i,1:5,j));
            
         psi.delta.sum_radians(i,1,j) = (psi.delta.sum(i,1,j)/180)*pi;
            
         psi.delta.table(i,1,j) = (TableA4_3(i,6) + TableA4_3(i,7) .*...
              date.JCE(j)) .* sin(psi.delta.sum_radians(i,1,j));
            
         epsilon.delta.table(i,1,j) =(TableA4_3(i,8)+TableA4_3(i,9).*...
              date.JCE(j)) .* cos(psi.delta.sum_radians(i,1,j));
    end
end

%Total delta PSI (nutation in longitude) (degrees) NREL equation (22)

psi.delta = sum(psi.delta.table)/36000000;

%Total delta EPSILON (nutation in obliquity) (degrees) NREL equation (23)

epsilon.delta = sum(epsilon.delta.table)./36000000;

%% CALCULATE TRUE OBLIQUITY OF THE ECLIPTIC (EPSILON) (DEGREES)

%Calculate mean obliquity of the ecliptic (EPSILON naught) 
%(in arc seconds) NREL equation (24)

U = (date.JME)/10;

epsilon.naught = 84381.448 + U .* (-4680.93 + U .*(-1.55 + U .* ...
    (1999.25 + U .*(-51.38 + U .*(-249.67 + U .* (-39.05 + U .* ...
    (7.12 + U .*(27.87 + U .* (5.79 + U .* 2.45)))))))));

%True obliquity of the ecliptic (EPSILON) (degrees) NREL equation (25)

for i = 1:length(epsilon.naught)
    for k = 1:length(date.JCE)
        epsilon.total(i,1) = ((epsilon.naught(i)./3600)) + ...
            epsilon.delta(:,:,k);
    end
end

%% CALCULATE THE ABERATTION CORRECTION (DELTA TAU)

%Delta TAU (degrees) NREL equation (26)

tau.delta = -20.4898./(3600.*R.AU);

%% CALCULATE THE APPARENT SUN LONGITUDE (LAMBDA)

%LAMBDA (degrees) NREL equation (27)

for i = 1:length(geocentric.normalized.longitude)
    for k = 1:length(date.JME)
        lambda.total(i) = geocentric.normalized.longitude(i) + ...
            psi.delta(:,:,k) + tau.delta(i);
    end
end


%% CALCULATE THE APPARENT SIDEREAL TIME AT GREENWICH AT ANY TIME GIVEN (NU)

%Mean sidereal time at Greenwich, (NU NAUGHT) (degrees) NREL equation (28)

nu.naught = 280.46061837 + 360.98564736629 .* ((date.jd)-2451545) ...
    + 0.000387933 .* (date.JC).^2 - (((date.JC).^3) ./ 38710000);

%Limite NU naught as described in 3.2.6

limit.M = abs(rem(nu.naught/360,1));

for i = 1:length(nu.naught)
    if nu.naught(i) >= 0
        nu.normalized(i) = 360*limit.M(i);
    elseif nu.naught(i) < 0 
        nu.normalized(i) = 360-360*limit.M(i);
    end
end

%Apparent sidereal time at Greenwich (NU) (degrees) NREL equation (29)

for i = 1:length(nu.normalized)
    for k = 1:length(date.JME)
        epsilon.rad(i) = (epsilon.total(i)/180)*pi;
        
        nu.total(i) = (nu.normalized(i)) + (psi.delta(:,:,k)) .* ...
            (cos(epsilon.rad(i)));
    end
end


%% CALCULATE THE GEOCENTRIC SUN RIGHT ASCENSION (ALPHA) 

%Sun right ascension (alpha) (radians) NREL equation (30)

for i = 1:length(lambda.total)
    lambda.rad(i) = (lambda.total(i)/180)*pi;
    geocentric.latitude_rad(i) = (geocentric.latitude(i)/180)*pi;
    
    alpha.radians(i,1) = atan2(sin(lambda.rad(i)) .* ...
        cos(epsilon.rad(i)) - tan(geocentric.latitude_rad(i)) .* ...
        sin(epsilon.rad(i)), cos(lambda.rad(i)));
end

%ALPHA in degrees using NREL section 3.2.5

alpha.degrees = (alpha.radians*180)/pi;

%ALPHA normalized to 360 degrees using NREL section 3.2.6

limit.Z = abs(rem(alpha.degrees/360,1));

for i = 1:length(alpha.degrees)
    if alpha.degrees(i) >= 0
        alpha.normalized(i) = 360*limit.Z(i);
    elseif alpha.degrees(i) < 0 
        alpha.normalized(i) = 360-360*limit.Z(i);
    end
end
  
%% CALCULATE THE GEOCENTRIC SUN DECLINATION (DELTA) 

%Sun declination (DELTA) (radians) NREL equation (31)

for i = 1:length(geocentric.latitude)
    
    delta.radians(i,1) = (asin(sin(geocentric.latitude_rad(i)) .* ...
        cos(epsilon.rad(i)) + cos(geocentric.latitude_rad(i)) .* ...
        sin(epsilon.rad(i)) .* sin(lambda.rad(i))));
end

%Convert to degrees

delta.degrees = (delta.radians*180)/pi;

%% CALCULATE THE OBSERVER LOCAL HOUR ANGLE (CAPITAL ETA, H) 

%Observer local hour angle (ETA) (degrees) NREL equation (32)

eta.degrees = nu.total + Inputs.longitude - alpha.normalized;

%ETA normalized to 360 degrees using NREL section 3.2.6

limit.Y = abs(rem(eta.degrees/360,1));

for i = 1:length(eta.degrees)
    if eta.degrees(i) >= 0
        eta.normalized(i) = 360*limit.Y(i);
    elseif eta.degrees(i) < 0 
        eta.normalized(i) = 360-360*limit.Y(i);
    end
end  

%% CALCULATE THE TOPOCENTRIC SUN RIGHT ASCENSION (ALPHA PRIME) 

%Equatorial horizontal parallax of the sun (XI) (degrees) 
%NREL equation (33)

xi.total = 8.794 ./ (3600.*R.AU);

%J (radians) NREL equation (34)

u.total = atan(0.99664719 .* tan((Inputs.latitude./180) .* pi));

%x NREL equation (35) (radians)

x.total = cos(u.total) + (((Inputs.elevation) ./ 6378140) .* ...
    cos((Inputs.latitude ./ 180).*pi));

%y NREL equation (36) (radians)

y.total = 0.99664719 .* (sin(u.total)) + ...
    (((Inputs.elevation) ./ 6378140) .* ...
    sin((Inputs.latitude ./ 180) .* pi));

%Parallax in the sun right ascension (DELTA ALPHA) (degrees) 
%NREL equation(37)

alpha.delta.radians = atan2(((-x.total) .* ...
    (sin((xi.total ./ 180) .* pi)) .* ...
    (sin((eta.normalized ./ 180) .*pi))), ...
    ((cos(delta.radians')) - (x.total) .* (sin((xi.total ./ 180) .*pi)) ...
    .* (cos((eta.normalized ./ 180) .* pi))));

%Convert to degrees

alpha.delta.degrees = (alpha.delta.radians*180)/pi;

%Topocentric sun right ascension (alpha prime) (degrees) NREL equation (38)

alpha.prime = alpha.normalized + alpha.delta.degrees;

%Topocentric sun declination (delta prime) (radians) NREL equation (39)

delta.prime.radians = atan2(((sin(delta.radians') - (y.total) .* ...
    (sin((xi.total ./ 180) .* pi))) .* (cos(alpha.delta.radians))), ...
    ((cos(delta.radians')) - ((x.total) .* (sin((xi.total./180).*pi)).*...
    (cos((eta.normalized ./ 180) .* pi)))));

%Convert to degrees

delta.prime.degrees = (delta.prime.radians .* 180) ./ pi;

%% CALCULATE THE TOPOCENTRIC LOCAL HOUR ANGLE (ETA PRIME) 

%Topocentric local hour angle (degrees) NREL equation (40)

eta.prime = eta.normalized - alpha.delta.degrees;

%% CALCULATE THE TOPOCENTRIC ZENITH ANGLE (THETA) 

%Topocentric elevation angle without atmospheric refraction correction 
%(e naught) (radians) NREL equation (41)

e.naught.radians = asin((sin((Inputs.latitude ./ 180) .* pi)) .* ...
    (sin(delta.prime.radians)) + (cos((Inputs.latitude ./ 180) .*pi)).* ...
    (cos(delta.prime.radians)) .* cos((eta.prime ./ 180) .* pi));

%Convert to degrees

e.naught.degrees = (e.naught.radians .* 180) ./ pi;

%Atmospheric refraction correction (delta e) (degrees) NREL equation (42)

e.delta.degrees = (Inputs.pressure ./ 1010) .* ...
    (283 ./ (273 + Inputs.temperature)) .* (1.02 ./ (60 .* ...
    tan(e.naught.radians + 10.3 ./ (e.naught.radians + 5.11))));

%Topocentric elevation angle (e) (degrees) NREL equation (43)

e.degrees = e.naught.degrees + e.delta.degrees;

%Topocentric zenith angle (theta) (degrees) NREL equation (44)

theta.degrees = 90 - e.degrees;

%% CALCULATE THE TOPOCENTRIC AZIMUTH ANGLE (PHI) 

%Topocentric astronomers azimuth angle (gamma) (radians) NREL equation (45)

gamma.radians = atan2((sin((eta.prime ./ 180) .* pi)), ...
    (((cos((eta.prime ./ 180) .* pi)) .* ...
    (sin((Inputs.latitude ./ 180) .* pi))) ...
    - ((tan(delta.prime.radians)) .* (cos((Inputs.latitude ./180).*pi)))));

%Convert to degrees

gamma.degrees = (gamma.radians*180)/pi;

%GAMMA normalized to 360 degrees using NREL section 3.2.6

limit.W = abs(rem(gamma.degrees/360,1));

for i = 1:length(gamma.degrees)
    if gamma.degrees(i) >= 0
        gamma.normalized(i) = 360*limit.W(i);
    elseif gamma.degrees(i) < 0 
        gamma.normalized(i) = 360-360*limit.W(i);
    end
end  

%Topocentric azimuth angle (phi) (degrees)

phi.degrees = gamma.normalized + 180;

%PHI normalized to 360 degrees using NREL section 3.2.6

limit.S = abs(rem(phi.degrees/360,1));

for i = 1:length(phi.degrees)
    if phi.degrees(i) >= 0
        phi.normalized(i) = 360*limit.S(i);
    elseif phi.degrees(i) < 0 
        phi.normalized(i) = 360-360*limit.S(i);
    end
end  

%% Output variables

%Set output variables to necessary values
Solar_Azimuth = phi.normalized;
Zenith = theta.degrees;

%% Calculate Incidence angle for SAM

%Local Standard Meridian
LSTM = (ceil(15.*(Inputs.longitude./15)));

%Angle of incidence between direct irradiance and solar panel
Incidence_Angle = acosd(cosd(90.*ones(size(LSTM))-Zenith) ...
    .*cosd(Solar_Azimuth-Panel_Azimuth).* ...
    sind(Panel_Tilt_Bound)+sind(90.*ones(size(LSTM))-Zenith).* ...
    cosd(Panel_Tilt_Bound)); 

%ensures no negative cos(AoI) values; caps loss to "negative" collection values
Incidence_Angle(Incidence_Angle > 90 & Incidence_Angle < 270) = 90;  
end



