function [Clear_Sky_Direct_Spectrum, Total_Clear_Sky_Direct, Clear_Sky_GHI_Spectrum, Total_Clear_Sky_GHI ] = Solar_Spectra(Latitude, Longitude, Day, Zenith)
% This function calculates the clear sky solar spectrum and total
% irradiance for a given lat, long, and time

%% Solar spectrum (Bird Model)

Lambda = [0.300 0.305 0.310 0.315 0.320 0.325 0.330 0.335 0.340 0.345 0.350 0.360 0.370 0.380 0.390 0.400 0.410 0.420 0.430 0.440...
    0.450 0.460 0.470 0.480 0.490 0.500 0.510 0.520 0.530 0.540 0.550 0.570 0.593 0.610 0.630 0.656 0.668 0.690 0.710 0.718...
    0.724 0.740 0.753 0.758 0.763 0.768 0.780 0.800 0.816 0.824 0.832 0.840 0.860 0.880 0.905 0.915 0.925 0.930 0.937 0.948...
    0.965 0.980 0.994 1.040 1.070 1.100 1.120 1.130 1.145 1.161 1.170 1.200 1.240 1.270 1.290 1.320 1.350 1.395 1.443 1.463...
    1.477 1.497 1.520 1.539 1.558 1.578 1.592 1.610 1.630 1.646 1.678 1.740 1.800 1.860 1.920 1.960 1.985 2.005 2.035 2.065...
    2.100 2.148 2.198 2.270 2.360 2.450 2.500 2.600 2.700 2.800 2.900 3.000 3.100 3.200 3.300 3.400 3.500 3.600 3.700 3.800...
    3.900 4.000]; %Wavelengths [micro(mu) m]

alpha = 1.2060.*ones(size(Lambda)); %alpha depends on the value of Lambda
alpha(Lambda <= 0.5) = 1.0274;  %larger wavelengths alpha = 1.2060 and smaller wavelengths alpha = 1.0274

Cs = ones(size(Lambda));    % Scattered irradiance correction factor [-]
Cs(Lambda < 0.45) = (Lambda(Lambda < 0.45)+0.55).^(1.8);    % Scattered irradiance correction factor
Cs = repmat(Cs',1,length(Zenith));

AM0 = [535.9 558.3 622 692.7 715.1 832.9 961.9 931.9 900.6 911.3 975.5 975.9 1119.9 1103.8 1033.8 1479.1 1701.3 1740.4 1587.2 1837 ...
    2005 2043 1987 2027 1896 1909 1927 1831 1891 1898 1892 1840 1768 1728 1658 1524 1531 1420 1399 1374 1373 1298 1269 1245 1223 ...
    1205 1183 1148 1091 1062 1038 1022 998.7 947.2 893.2 868.2 829.7 830.3 814 786.9 768.3 767 757.6 688.1 640.7 606.2 585.9 570.2 ...
    564.1 544.2 533.4 501.6 477.5 442.7 440 416.8 391.4 358.9 327.5 317.5 307.3 300.4 292.8 275.5 272.1 259.3 246.9 244 243.5 234.8 ...
    220.5 190.8 171.1 144.5 135.7 123 123.8 113 108.5 97.5 92.4 82.4 74.6 68.3 63.8 49.5 48.5 38.6 36.6 32 28.1 24.8 22.1 19.6 17.5 ...
    15.7 14.1 12.7 11.5 10.4 9.5  8.6]; AM0 = repmat(AM0',1,length(Zenith)); %Extraterrestrial Irradiance at each wavelength (Lambda)[W/m^2/?m]

aw = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.075 0 0 0 0 0.016 0.0125 1.8 2.5 0.061 0.0008 0.0001 0.00001 0.00001...
    0.0006 0.036 1.6 2.5 0.5 0.155 0.00001 0.0026 7 5 5 27 55 45 4 1.48 0.1 0.00001 0.001 3.2 115 70 75 10 5 2 0.002 0.002 0.1 4 200 1000 ...
    185 80 80 12 0.16 0.002 0.0005 0.0001 0.00001 0.0001 0.001 0.01 0.036 1.1 130 1000 500 100 4 2.9 1 0.4 0.22 0.25 0.33 0.5 4 80 310 15000 ...
    22000 8000 650 240 230 100 120 19.5 3.6 3.1 2.5 1.4 0.17 0.0045]; aw = repmat(aw',1,length(Zenith)); %Absorption spectrum of water vapor

ao = [10 4.8 2.7 1.35 0.8 0.38 0.16 0.075 0.04 0.019 0.007 0 0 0 0 0 0 0 0 0 0.003 0.006 0.009 0.014 0.021 0.03 0.04 0.048 0.063 0.075 0.085...
    0.12 0.119 0.12 0.09 0.065 0.051 0.028 0.018 0.015 0.012 0.01 0.008 0.007 0.006 0.005 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Absorption spectrum of ozone

au = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.15 0 0 0 0 0 0 4 0.35 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0.05 0.3 0.02 0.0002 0.00011 0.00001 0.05 0.011 0.005 0.0006 0 0.005 0.13 0.04 0.06 0.13 0.001 0.0014 0.0001 0.00001 0.00001 ...
    0.0001 0.001 4.3 0.2 21 0.13 1 0.08 0.001 0.00038 0.001 0.0005 0.00015 0.00014 0.00066 100 150 0.13 0.0095 0.001 0.8 1.9 1.3 0.075 0.01 ...
    0.00195 0.004 0.29 0.025]; au = repmat(au',1,length(Zenith)); %Absorption spectrum of mixed gases

Po = 1013; %Corrective pressure in milibar
P = 1013; %measured surface pressure in milibar
R = P/Po; %Pressure ratio
beta = .125; % Constant for turbidity calculations
W = 4; %precipitable water vapor in cm
ho = 22000; %height of Ozone concentration
ALG = log(1-0.65);
AFS = ALG*(1.459+ALG*(0.1595+ALG*0.4129));
BFS = ALG*(0.0783+ALG*(-0.3824-ALG*0.5874));
Fs_Prime = 1-0.5*exp((AFS+BFS/1.8)/1.8);  % Fraction of aerosol scatter
rg = 0.2; %Ground Albedo

DA = 2*pi*((Day-1)/365); %day angle in radians
D = 1.00011+0.034221*cos(DA)+0.00128*sin(DA)+0.000719*cos(2*DA)+0.000077*sin(2*DA); %earth-sun distance factor

o3 = (235+(150+40*sin(0.9856*(Day-30))+20.*sin(3.*(Longitude+20))).*(sin(1.28.*Latitude).^2))./1000; %Amount of Ozone in atmosphere centimetres
M =(cosd(Zenith)+0.15.*(93.885-Zenith).^(-1.253)).^-1; M = repmat(M',length(Lambda),1);  %relative air mass
M_Prime = M.*R; %pressure-corrected air mass
Mo = (1+(ho/6370))./(((cos(Zenith)).^2+(2*ho/6370)).^(0.5)); %Ozone mass
Fs = 1 - 0.5.*exp((AFS+BFS.*cosd(Zenith)).*cosd(Zenith));  % Fraction of aerosol scatter dependent of zenith angle

% Direct Irradiance Calculations
Tr = exp(-M_Prime./((repmat(Lambda',1,length(Zenith)).^4).*(115.6406-(1.335./(repmat(Lambda',1,length(Zenith)).^2))))); % Rayleigh Scattering
Tr_Prime = Tr.*R; % Corretced Rayleigh Scattering
TAUa = (beta).*(Lambda.^(-alpha)); TAUa = repmat(TAUa',1,length(Zenith));  % Aerosol extinction
Ta = exp(-TAUa.*M); % Aerosol Scattering and Absorption
Tw = exp((-0.2385.*aw.*W.*M)./(1+20.07.*aw.*W.*M).^(0.45)); % Water Vapor Absorption
Tw_Prime = Tw.*R; % Corrected Water Vapor Absorption
To = exp(-repmat(ao',1,length(Zenith)).*repmat(o3',length(Lambda),1).*repmat(Mo',length(Lambda),1)); % Ozone Absorption
To_Prime = To.*R; % Corrected Ozone Absorption
Tu = exp((-1.41.*au.*M_Prime)./((1+118.93.*au.*M_Prime).^(0.45))); %Uniformly mixed gas absorption
Clear_Sky_Direct_Spectrum = AM0.*D.*Tr.*Ta.*Tw.*To.*Tu; %Clear sky direct irradiance spectrum

% Diffuse Irradiance Calculations
w = 0.945.*exp(-0.095.*(log(Lambda./0.4)).^2); w = repmat(w',1,length(Zenith));
Taa = exp(-(1-w).*TAUa.*M); % Transmittance of aerosol absorption
Taa_Prime = Taa.*R; % Corrected transmittance of aerosol absorption
Tas = exp(-w.*TAUa.*M); % Transmittance of aerosol scattering
Tas_Prime = Tas.*R; % Corrected transmittance of aerosol scattering
rs = (To_Prime.*Tw_Prime.*Taa_Prime).*(0.5.*(1-Tr_Prime)+(1-Fs_Prime).*Tr_Prime.*(1-Tas_Prime)); % Ground albedo function
Ir = AM0.*D.*repmat(cosd(Zenith'),length(Lambda),1).*To.*Tu.*Tw.*Taa.*(1-Tr.^0.95).*0.5.*Cs; % Irradiance through Rayleigh scattering
Ia = AM0.*D.*repmat(cosd(Zenith'),length(Lambda),1).*To.*Tu.*Tw.*Taa.*(Tr.^1.5).*(1-Tas).*repmat(Fs',length(Lambda),1).*Cs; % Irradiance through aerosol scattering
Ig = (Clear_Sky_Direct_Spectrum.*repmat(cosd(Zenith'),length(Lambda),1)+Ir+Ia).*rs.*rg.*Cs./(1-rs.*rg); % Irradiance from reflection of ground and air
Clear_Sky_Diffuse_Spectrum = Ir+Ia+Ig; %Clear sky diffuse irradiance spectrum

Clear_Sky_GHI_Spectrum = Clear_Sky_Direct_Spectrum.*repmat(cosd(Zenith'),length(Lambda),1) + Clear_Sky_Diffuse_Spectrum; % Clear sky GHI spectrum [W/m^2/?m]

%% Direct and Diffuse Total

Total_Clear_Sky_Direct = trapz(Lambda(16:end),Clear_Sky_Direct_Spectrum(16:end,:)); % Total direct irradiance within range CERES instruments can sense [W/m^2]
Total_Clear_Sky_GHI = trapz(Lambda(16:end),Clear_Sky_GHI_Spectrum(16:end,:)); % Total GHI within range CERES instruments can sense [W/m^2]

end

