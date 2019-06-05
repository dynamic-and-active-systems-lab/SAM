function [Zenith,Incidence_Angle] = Solar_Position(Latitude, Longitude, Day_of_Year, Minute, Panel_Tilt_Bound, Panel_Azimuth)
% This function calculates the clear sky solar spectrum and total
% irradiance for a given lat, long, and time

%% Solar Position Calculator

LSTM = (ceil(15.*(Longitude./15))); % Local Standard Meridian
B = (360/365)*(Day_of_Year+81).*ones(size(LSTM)); % Coefficient for Equation of Time
ET = 9.87.*sind(2.*B)-7.53.*cosd(B)-1.5.*sind(B); % Equation of Time in minutes
Solar_Declination = 23.45.*sind(((360*(284+Day_of_Year))/365)).*ones(size(LSTM)); % Solar Declination

% %get rid of nans i Longitude without changing size of matrix; THIS IS BAD
% Longitude(isnan(Longitude)) = max(max(Longitude))-90; %set it so time zone is very off from real measurements
% 
% TZ_offset = (-timezone(Longitude)-Time_Zone); %time zone offset
% TZ_offset = reshape(TZ_offset, 86400, []);
% 
% a = -timezone(max(max(Longitude))+90)-Time_Zone; %find our impossible value locations
% TZ_offset(TZ_offset == a) = nan; %set back to nan

Solar_Time = Minute + (4.*(LSTM - Longitude)) + ET; % Solar Time in minutes
Hour_Angle = (Solar_Time-720.*ones(size(LSTM)))./4; % Hour angle in minutes
Zenith = 90.*ones(size(LSTM))-(asind(sind(Latitude).*sind(Solar_Declination)+cosd(Latitude).*cosd(Solar_Declination).*cosd(Hour_Angle))); % Calculates solar zenith angle
Zenith(Zenith > 90) = 90;
Solar_Azimuth = sign(Hour_Angle).*acosd((cosd(Zenith).*sind(Latitude)-sind(Solar_Declination))./(sind(Zenith).*cosd(Latitude)))+180.*ones(size(LSTM)); % Solar azimuth angle in degrees north
Incidence_Angle = acosd(cosd(90.*ones(size(LSTM))-Zenith).*cosd(Solar_Azimuth-Panel_Azimuth).*sind(Panel_Tilt_Bound)+sind(90.*ones(size(LSTM))-Zenith).*cosd(Panel_Tilt_Bound)); % Angle of incidence between direct irradiance and solar panel
Incidence_Angle(Incidence_Angle > 90 & Incidence_Angle < 270) = 90;   %ensures no negative cos(AoI) values; caps loss to "negative" collection values

end