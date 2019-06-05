function [Total_All_Sky_Direct_Down, Total_All_Sky_GHI] = Irradiance_ReadIn_Week(Year, Month, Start, End, Week_Length, Time_Zone, Day_of_Month, Latitude_CERES_Indices, Longitude_CERES_Indices, irr_file_path)
%% Functions Calculations
Month_str = {'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'};  %used for irradiance variable name read in
Month_Length = [31 28 31 30 31 30 31 31 30 31 30 31;...
    31 29 31 30 31 30 31 31 30 31 30 31]; %# of days in month, row 1: Common Year row 2: Leap Year
Month_Start_Day = [1 32 60 91 121 152 182 213 244 274 305 335;...
    1 32 61 92 122 153 183 214 245 275 306 336]; %start day of each month, row 1: Common Year row 2: Leap Year
Month_End_Day = Month_Start_Day + Month_Length - 1;
row = (mod(Year,4) == 0)*2 + (mod(Year,4) ~= 0)*1;  %logically determine row
% Day_of_Month = (Start-Month_Start_Day(row,Month)) + 1;    %find day of the month of current day

%check if Time Zone input was valid
if abs(Time_Zone) > 12 %outside range of possible time zones
    error('Invalid time zone, please check ''Time Zone'' input and try again.');    %throw error and stop model
end

if Time_Zone < 0 %negative time zone (UTC - __)
    if End == Month_End_Day(row,Month) %if month deployment reaches last day of the month
        %set variable for next month & year; accounting for possible change of calendar year
        if Month == 12 %if current month is Dec
            nextMonth = 1; nextYear = Year + 1; %next month is Jan and iterate year by one
        else %if current month is anything but Dec
            nextMonth = Month + 1; nextYear = Year; %iterate month by one and keep year the same
        end       
        %read in current and next month of irradiance data
        org_file_path = cd(irr_file_path); %change directory according to file holding Irradiance data
        load(sprintf('CERES_SYN1deg_1H_TrAq_Ed4a_%d_%s_SW_Irr_Flux.mat',Year,string(Month_str(Month,:))))   %load current month's irradiance data
        TAS_Dir_h = Total_All_Sky_Direct_h; %create dummy variable with this months irradiance data for interpolation
        TAS_GHI_h = Total_All_Sky_GHI_h;
        load(sprintf('CERES_SYN1deg_1H_TrAq_Ed4a_%d_%s_SW_Irr_Flux.mat',nextYear,string(Month_str(nextMonth,:))))   %load next month's irradiance data
        TAS_Dir_h = cat(3,TAS_Dir_h,Total_All_Sky_Direct_h);    %append variable with next months irradiance data
        TAS_GHI_h = cat(3,TAS_GHI_h,Total_All_Sky_GHI_h);
        clear Total_All_Sky_Direct_h Total_All_Sky_Diffuse_h
        cd(org_file_path) %change to original directory  
    else %if deployment does not reach last day of month
        %read in current month irradiance data
        org_file_path = cd(irr_file_path); %change directory according to file holding Irradiance data
        load(sprintf('CERES_SYN1deg_1H_TrAq_Ed4a_%d_%s_SW_Irr_Flux.mat',Year,string(Month_str(Month,:))))   %load current month's irradiance data
        TAS_Dir_h = Total_All_Sky_Direct_h; %create dummy variable with this months irradiance data, for interpolation
        TAS_GHI_h = Total_All_Sky_GHI_h;
        cd(org_file_path) %change to original directory
    end
    %create start/end indices corresponding to ceres irradiance matrix time indices
    ceres_strt = (Day_of_Month - 1)*24 + 1; %first ceres irr time index of month deployment in UTC
    ceres_end = ceres_strt + (Week_Length*24) - 1; %last ceres irr time index of current month in UTC
    
elseif Time_Zone > 0 %positive time zone (UTC + __)
    if sum(Start == Month_Start_Day(row,:)) == 1 %if the first day of the month is experienced
        %set variable for previous month & year; accounting for possible change of calendar year
        if Month == 1 %if it is Jan
            prevMonth = 12; prevYear = Year - 1; %prev month is Dec and decrease year by one
        else %if it is not Jan
            prevMonth = Month - 1; prevYear = Year; %decrease month by one and keep year the same
        end
        %read in current and previous month of irradiance data and use for interpolation
        org_file_path = cd(irr_file_path); %change directory according to file holding Irradiance data
        load(sprintf('CERES_SYN1deg_1H_TrAq_Ed4a_%d_%s_SW_Irr_Flux.mat',prevYear,string(Month_str(prevMonth,:))))   %load previous month's irradiance data
        prevMonth_len = size(Total_All_Sky_Direct_h,3); %get length of time indices for use in ceres_strt
        TAS_Dir_h = Total_All_Sky_Direct_h; %create dummy variable with this months irradiance data, for interpolation
        TAS_GHI_h = Total_All_Sky_GHI_h;
        clear Total_All_Sky_Direct_h Total_All_Sky_Diffuse_h
        load(sprintf('CERES_SYN1deg_1H_TrAq_Ed4a_%d_%s_SW_Irr_Flux.mat',Year,string(Month_str(Month,:))))   %load current month's irradiance data
        TAS_Dir_h = cat(3,TAS_Dir_h,Total_All_Sky_Direct_h);    %append variable with next months irradiance data
        TAS_GHI_h = cat(3,TAS_GHI_h,Total_All_Sky_GHI_h);
        clear Total_All_Sky_Direct_h Total_All_Sky_Diffuse_h
        cd(org_file_path) %change to original directory
        %create the first ceres irr time index of month deployment in UTC
        ceres_strt = prevMonth_len + 1; %when ceres start depends on length of previous month 
        
    else  %if it is not the first day of the month
        %read in current month of irradiance data
        org_file_path = cd(irr_file_path); %change directory according to file holding Irradiance data
        load(sprintf('CERES_SYN1deg_1H_TrAq_Ed4a_%d_%s_SW_Irr_Flux.mat',Year,string(Month_str(Month,:))))   %load current month's irradiance data
        TAS_Dir_h = Total_All_Sky_Direct_h; %create dummy variable with this months irradiance data, for interpolation
        TAS_GHI_h = Total_All_Sky_GHI_h;
        cd(org_file_path) %change to original directory
        %create the first ceres irr time index of month deployment in UTC
        ceres_strt = (Day_of_Month - 1)*24 + 1; %when ceres start does NOT depend on length of previous month 
        
    end
    %create vector corresponding to ceres irradiance matrix time indices
    ceres_end = ceres_strt + (Week_Length*24) - 1; %last ceres irr time index of current month in UTC
    
else %zero time zone (UTC)
    org_file_path = cd(irr_file_path); %change directory according to file holding Irradiance data
    %read in current month of irradiance data
    load(sprintf('CERES_SYN1deg_1H_TrAq_Ed4a_%d_%s_SW_Irr_Flux.mat',Year,string(Month_str(Month,:))))   %load current month's irradiance data
    TAS_Dir_h = Total_All_Sky_Direct_h; %create dummy variable with this months irradiance data, for interpolation
    TAS_GHI_h = Total_All_Sky_GHI_h;
    cd(org_file_path) %change to original directory
    
    %create vector corresponding to ceres irradiance matrix time indices
    ceres_strt = (Day_of_Month - 1)*24 + 1; %first ceres irr time index of month deployment in UTC
    ceres_end = ceres_strt + (Week_Length*24) - 1; %last ceres irr time index of current month in UTC

end
%%
%%%%%% Interpolate
%create vector corresponding to ceres irradiance matrix time indices accounting for time zone
ceres_ind = (ceres_strt:ceres_end)-Time_Zone; 
%create time vector for interpolation
t_ceres = (0.5:(Week_Length*24 - 0.5))*3600; %time in ceres irradiance matrix resolution (0.5, 1.5, ..., 32.5 for every day of month), but in units of seconds
t_seconds = 1:86400*Week_Length; %time of deployment in seconds units and second resolution
%find lon and lat min and max of month
Lat_wk_max = max(max(Latitude_CERES_Indices));
Lat_wk_min = min(min(Latitude_CERES_Indices));
Lon_wk_max = max(max(Longitude_CERES_Indices));
Lon_wk_min = min(min(Longitude_CERES_Indices));
%create new indices to use in main code
% Latitude_Month_Indices = Latitude_CERES_Indices - Lat_min + 1;
% Longitude_Month_Indices = Longitude_CERES_Indices - Lon_min + 1;
%temp variables for pulling correct indices in interpolation
A = Lon_wk_min:Lon_wk_max;    
B = Lat_wk_min:Lat_wk_max;
TAS_Dir_interp_vect = permute(TAS_Dir_h(A, B, ceres_ind), [3 1 2]);
TAS_GHI_interp_vect = permute(TAS_GHI_h(A, B, ceres_ind),[3 1 2]);
a = round(linspace(1, length(Lon_wk_min:Lon_wk_max+1), 4));
b = round(linspace(1, length(Lat_wk_min:Lat_wk_max+1), 4));
%preallocate matrix (seconds of deployment in month x degrees of lon in month x degrees of lat in month)
Total_All_Sky_Direct_Down = zeros(length(t_seconds), length(A), length(B)); 
Total_All_Sky_GHI = zeros(length(t_seconds), length(A), length(B));
if length(A) > 10 && length(B) > 10 %if monthly lat and lon are large
    for i = 1:3
        for j = 1:3
            Total_All_Sky_Direct_Down(:, a(i):a(i+1)-1, b(j):b(j+1)-1) = interp1(t_ceres, TAS_Dir_interp_vect(:, a(i):a(i+1)-1, b(j):b(j+1)-1), t_seconds);
            Total_All_Sky_GHI(:, a(i):a(i+1)-1, b(j):b(j+1)-1) = interp1(t_ceres, TAS_GHI_interp_vect(:, a(i):a(i+1)-1, b(j):b(j+1)-1), t_seconds);
        end
    end
else
    Total_All_Sky_Direct_Down = interp1(t_ceres, TAS_Dir_interp_vect, t_seconds);
    Total_All_Sky_GHI = interp1(t_ceres, TAS_GHI_interp_vect, t_seconds);
end
    

%%%%%% Reshape Matrix Before Passing Back to MAIN Function
Total_All_Sky_Direct_Down = permute(Total_All_Sky_Direct_Down, [2 3 1]);
Total_All_Sky_GHI = permute(Total_All_Sky_GHI, [2 3 1]);
Total_All_Sky_Direct_Down = reshape(Total_All_Sky_Direct_Down,length(A),length(B),86400,[]); %dimensions of relative month lon x relative month lat x seconds x day
Total_All_Sky_GHI = reshape(Total_All_Sky_GHI,length(A),length(B),86400,[]);

end