 function [out, wb_progress] = MAIN(wb, wb_progress, wb_str, wb_unit, filename, Water_Quality, Time_Zone, outputNames, in_file_path, irr_file_path, array_type, array_n, Solar_Cell_Area, Rsh, Rs, C, Eg_0, alpha, beta)
%% Load in data
waitbar(wb_progress,wb,[wb_str,'Loading Data...']);

load(sprintf('%s/%s', in_file_path, filename)); % Reads Start and End Day

%% Time of year
waitbar(wb_progress,wb,[wb_str,'Initializing Calculations...']);

Month_Start_Day = [1 32 60 91 121 152 182 213 244 274 305 335;...
    1 32 61 92 122 153 183 214 245 275 306 336]; %start day of each month, row 1: Common Year row 2: Leap Year
 
dateTime = datetime(Time,'ConvertFrom','excel');
Year = unique(dateTime.Year,'stable');
Month = unique(dateTime.Month,'stable');

row = (mod(Year,4) == 0)*2 + (mod(Year,4) ~= 0)*1;  %logically determine row of 'Month_Start_Day' vector
Day_of_Year =  unique(dateTime.Day,'stable') + Month_Start_Day(row, Month) - 1;
Start = Day_of_Year(1); % Start day of data
End = Day_of_Year(end); % End day of data
File_Deployment_Length = End-Start+1;

%% Tag Variables

Depth = reshape(Depth,86400,File_Deployment_Length); % Reshapes matrix to 86400 by # days of data
Depth(Depth < 0) = 0;   %removes negative enteries

T = reshape(T,86400,File_Deployment_Length) + 273.15; %temperature of cell in K

Panel_Azimuth = reshape(Panel_Azimuth,86400,File_Deployment_Length); % Azimuth angle of panel in degrees north (values 0-359)
Panel_Tilt = reshape(Panel_Tilt,86400,File_Deployment_Length); % Tilt of solar panel (values 0-90)
Panel_Tilt_Bound = Panel_Tilt;
Panel_Tilt_Bound(Panel_Tilt > 90 & Panel_Tilt < 270) = 90; %PT in this range will cause negative power, when it should just cause zero power output

%% Latitude and Longitude

Latitude = reshape(Latitude,86400,File_Deployment_Length); % Negative for South and postive for North

Longitude = reshape(Longitude,86400,File_Deployment_Length); % Negative for West and postive for East

%convert Latitude and Longitude values into indices corresponding to
%indices used by CERES Irradiance data
isneg = @(val) val < 0; %set up quick functions for clean conversions
ispos = @(val) val >= 0;
Latitude_CERES_Indices = round(Latitude) + 91 - (round(Latitude)== 90);
Longitude_CERES_Indices = round(Longitude) + ispos(round(Longitude)) + isneg(round(Longitude)).*361;

%find lat/lon indices relative to travelled in the month; matrix indices
Lat_month_max = max(max(Latitude_CERES_Indices)); %find lon and lat min and max of month
Lat_month_min = min(min(Latitude_CERES_Indices));
Lon_month_max = max(max(Longitude_CERES_Indices));
Lon_month_min = min(min(Longitude_CERES_Indices));
Latitude_Month_Indices = Latitude_CERES_Indices - Lat_month_min + 1;
Longitude_Month_Indices = Longitude_CERES_Indices - Lon_month_min + 1;
Lon_month = Lon_month_min:Lon_month_max;    
Lat_month = Lat_month_min:Lat_month_max;

%% Matrix Preallocation

Zenith = zeros(86400,File_Deployment_Length); % Preallocated matrix
Incidence_Angle = zeros(86400,File_Deployment_Length); % Preallocated matrix
Total_All_Sky_Direct_Down = zeros(length(Lon_month), length(Lat_month), 86400, File_Deployment_Length); %Preallocated matrix
Total_All_Sky_GHI = zeros(length(Lon_month), length(Lat_month), 86400, File_Deployment_Length); %Preallocated matrix
Clear_Sky_Direct_Spectrum = zeros(122,86400,File_Deployment_Length); % Preallocated matrix
Total_Clear_Sky_Direct = zeros(86400,File_Deployment_Length); % Preallocated matrix
Clear_Sky_GHI_Spectrum = zeros(122,86400,File_Deployment_Length); % Preallocated matrix
Total_Clear_Sky_GHI = zeros(86400,File_Deployment_Length); % Preallocated matrix
CERES_Direct_Down = zeros(86400,File_Deployment_Length); % Preallocated matrix
CERES_GHI = zeros(86400,File_Deployment_Length); % Preallocated matrix
Scaled_Direct_Down_Irradiance_Spectrum = zeros(122,86400,File_Deployment_Length); % Preallocated matrix
Scaled_Direct_Irradiance_Spectrum = zeros(122,86400,File_Deployment_Length); % Preallocated matrix
Scaled_Diffuse_Irradiance_Spectrum = zeros(122,86400,File_Deployment_Length); % Preallocated matrix
Scaled_GHI_Spectrum = zeros(122,86400,File_Deployment_Length); % Preallocated matrix
Total_Scaled_GHI = zeros(86400,File_Deployment_Length); % Preallocated matrix
Scaled_Transmitted_GHI_Spectrum = zeros(122,86400,File_Deployment_Length); % Preallocated matrix
Total_Scaled_Transmitted_GHI = zeros(86400,File_Deployment_Length); % Preallocated matrix
Total_Solar_Irradiance_Incident_To_Cell = zeros(86400,File_Deployment_Length); % Preallocated matrix
Solar_Cell_Power_Output = zeros(86400,File_Deployment_Length); % Preallocated matrix
V_oc = zeros(86400,File_Deployment_Length); % Preallocated matrix
I_sc = zeros(86400,File_Deployment_Length); % Preallocated matrix

%% Power output calculations

% Time = reshape(Time, 86400, []);

wk = round(linspace(1, File_Deployment_Length+1, 5)); %for week chunks
for i = 1:4 %calculations done a week at a time 
    % Solar Position calculator based on "Principles of Sustatinable Energy" by Kreith
    Minute = repmat(((1:86400)'./60), 1, length(wk(i):wk(i+1)-1));
    [Zenith(:,wk(i):wk(i+1)-1),Incidence_Angle(:,wk(i):wk(i+1)-1)] = Solar_Position(Latitude(:,wk(i):wk(i+1)-1), Longitude(:,wk(i):wk(i+1)-1), Day_of_Year(wk(i):wk(i+1)-1), Minute, Panel_Tilt_Bound(:,wk(i):wk(i+1)-1), Panel_Azimuth(:,wk(i):wk(i+1)-1)); %both angles [deg]
%     [Zenith(:,wk(i):wk(i+1)-1), Solar_Azimuth(:,wk(i):wk(i+1)-1), Incidence_Angle(:,wk(i):wk(i+1)-1)] = Solar_Angles_PST(Panel_Tilt_Bound(:,wk(i):wk(i+1)-1), Panel_Azimuth(:,wk(i):wk(i+1)-1), Latitude(:,wk(i):wk(i+1)-1), Longitude(:,wk(i):wk(i+1)-1), Time(:,wk(i):wk(i+1)-1));
    
    %Irradiance Read In
    a = min(min(Longitude_Month_Indices(:,wk(i):wk(i+1)-1))):max(max(Longitude_Month_Indices(:,wk(i):wk(i+1)-1)));
    b = min(min(Latitude_Month_Indices(:,wk(i):wk(i+1)-1))):max(max(Latitude_Month_Indices(:,wk(i):wk(i+1)-1)));
    Week_Length = wk(i+1)-wk(i); %(wk(i+1)-1) - wk(i) + 1 
    Day_of_Month = wk(i);
    [Total_All_Sky_Direct_Down(a,b,:,wk(i):wk(i+1)-1), Total_All_Sky_GHI(a,b,:,wk(i):wk(i+1)-1)] = Irradiance_ReadIn(Year, Month, Start, End, Week_Length, Time_Zone, Day_of_Month, Latitude_CERES_Indices(:,wk(i):wk(i+1)-1), Longitude_CERES_Indices(:,wk(i):wk(i+1)-1), irr_file_path); %Read in CERES data and interpolate only latitude and longitude values seen to a second resolution
    
end

for day = 1:File_Deployment_Length % daily calculations
    waitbar(wb_progress,wb,[wb_str,sprintf('Calculating Week %d/%d',day,File_Deployment_Length)]);
    
    % Solar Spectra calculations
    [Clear_Sky_Direct_Spectrum(:,:,day), Total_Clear_Sky_Direct(:,day), Clear_Sky_GHI_Spectrum(:,:,day), Total_Clear_Sky_GHI(:,day)] = Solar_Spectra(Latitude(:,day), Longitude(:,day), day, Zenith(:,day)); % Clear sky GHI and underwater GHI from Bird model; spectrums [W/m^2/?m], totals [W/m^2] 

    for sec = 1:86400
        if isnan(Latitude_Month_Indices(sec,day)) || isnan(Longitude_Month_Indices(sec,day))
            CERES_Direct_Down(sec,day) = nan;
            CERES_GHI(sec,day) = nan;
        else
            CERES_Direct_Down(sec,day) = Total_All_Sky_Direct_Down(Longitude_Month_Indices(sec,day),Latitude_Month_Indices(sec,day),sec,day);
            CERES_GHI(sec,day) = Total_All_Sky_GHI(Longitude_Month_Indices(sec,day),Latitude_Month_Indices(sec,day),sec,day);
        end
    end
    
    Direct_Down_Irradiance_Scale_Factor = CERES_Direct_Down(:,day)./Total_Clear_Sky_Direct(:,day);   
    Direct_Down_Irradiance_Scale_Factor(Total_Clear_Sky_Direct(:,day) == 0) = 0; %when denominator is 0 (at night) will create unrepresentative nan
    GHI_Scale_Factor = CERES_GHI(:,day)./Total_Clear_Sky_GHI(:,day);
    GHI_Scale_Factor(Total_Clear_Sky_GHI(:,day) == 0) = 0; %when denominator is 0 (at night) will create unrepresentative nan
    
    Scaled_Direct_Down_Irradiance_Spectrum(:,:,day) =  repmat(Direct_Down_Irradiance_Scale_Factor',122,1).*Clear_Sky_Direct_Spectrum(:,:,day); % Corrected direct irradiance spectrum
    Scaled_GHI_Spectrum(:,:,day) =  repmat(GHI_Scale_Factor',122,1).*Clear_Sky_GHI_Spectrum(:,:,day); % Scaled GHI spectrum
    
    Scaled_Diffuse_Irradiance_Spectrum(:,:,day) = Scaled_GHI_Spectrum(:,:,day) - Scaled_Direct_Down_Irradiance_Spectrum(:,:,day); % Calculates diffuse irradiance spectrum from direct irradiance and GHI spectra   
    Scaled_Diffuse_Irradiance_Spectrum(:,:,day) = Scaled_Diffuse_Irradiance_Spectrum(:,:,day).*(Scaled_Diffuse_Irradiance_Spectrum(:,:,day) > 0); %direct down irradiance can be greater than GHI, diffuse cannot be negative cause nature dude
    Scaled_Direct_Irradiance_Spectrum(:,:,day) = Scaled_Direct_Down_Irradiance_Spectrum(:,:,day)./repmat(cosd(Zenith(:,day)'),122,1);
    
    [Total_Scaled_GHI(:,day), Scaled_Transmitted_GHI_Spectrum(:,:,day), Total_Scaled_Transmitted_GHI(:,day)] = Transmitted_GHI(Scaled_Direct_Irradiance_Spectrum(:,:,day),Scaled_Diffuse_Irradiance_Spectrum(:,:,day),Scaled_GHI_Spectrum(:,:,day),Zenith(:,day));
    
    %Power Output Calculations
    [Total_Solar_Irradiance_Incident_To_Cell(:,day), Solar_Cell_Power_Output(:,day), V_oc(:,day),I_sc(:,day)] = Sub_Surface_Solar_Power(Scaled_Direct_Irradiance_Spectrum(:,:,day),Scaled_Diffuse_Irradiance_Spectrum(:,:,day), Scaled_Transmitted_GHI_Spectrum(:,:,day), Water_Quality, Depth(:,day), Panel_Tilt_Bound(:,day), Incidence_Angle(:,day),Total_Scaled_GHI(:,day),T(:,day), array_type, array_n, Solar_Cell_Area, Rsh, Rs, C, Eg_0, alpha, beta); % Total GHI and solar cell power output at each minute
    
    wb_progress = wb_progress + wb_unit;
end
Solar_Cell_Power_Output(isnan(Solar_Cell_Power_Output)) = 0;

%% Create Output Structure

out.Time = Time;
out.Depth = Depth;
out.Panel_Tilt = Panel_Tilt;
out.Panel_Tilt_Bound = Panel_Tilt_Bound;
out.Incidence_Angle = Incidence_Angle;
out.Water_Quality = Water_Quality;

for ioutput = 1:length(outputNames)
    outStr = sprintf('out.%s',outputNames{ioutput});
    outStr = [outStr, sprintf('=%s;',outputNames{ioutput})];
    eval(outStr)
end

 end