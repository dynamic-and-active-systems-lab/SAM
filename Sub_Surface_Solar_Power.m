function [Total_Solar_Irradiance_Incident_To_Cell, Solar_Cell_Power_Output, V_oc, IL] = Sub_Surface_Solar_Power(Scaled_Direct_Irradiance_Spectrum,Scaled_Diffuse_Irradiance_Spectrum, Scaled_Transmitted_GHI_Spectrum, Water_Quality, Depth, Panel_Tilt_Bound, Incidence_Angle,Total_Scaled_GHI, T, array_type, array_n, Solar_Cell_Area, Rsh, Rs, C, Eg_0, alpha, beta)
%This function calculates the power output of a solar cell placed at any
%depth underwater. 

%% Constant Definitions

% Spectral values
Lambda = [0.300 0.305 0.310 0.315 0.320 0.325 0.330 0.335 0.340 0.345 0.350 0.360 0.370 0.380 0.390 0.400 0.410 0.420 0.430 0.440...
      0.450 0.460 0.470 0.480 0.490 0.500 0.510 0.520 0.530 0.540 0.550 0.570 0.593 0.610 0.630 0.656 0.668 0.690 0.710 0.718...
      0.724 0.740 0.753 0.758 0.763 0.768 0.780 0.800 0.816 0.824 0.832 0.840 0.860 0.880 0.905 0.915 0.925 0.930 0.937 0.948...
      0.965 0.980 0.994 1.040 1.070 1.100 1.120 1.130 1.145 1.161 1.170 1.200 1.240 1.270 1.290 1.320 1.350 1.395 1.443 1.463...
      1.477 1.497 1.520 1.539 1.558 1.578 1.592 1.610 1.630 1.646 1.678 1.740 1.800 1.860 1.920 1.960 1.985 2.005 2.035 2.065...
      2.100 2.148 2.198 2.270 2.360 2.450 2.500 2.600 2.700 2.800 2.900 3.000 3.100 3.200 3.300 3.400 3.500 3.600 3.700 3.800...
      3.900 4.000]; %Wavelengths (micrometer)

EQE = [0 0.6447 0.6591 0.6731 0.6867 0.6888 0.6916 0.6940 0.6976 0.6982 0.6994 0.7018 0.7366 0.7734 0.8099 0.8435 0.86014 0.8760...
     0.8864 0.8982 0.9052 0.9124 0.9184 0.9258 0.9289 0.9314 0.9336 0.9367 0.9377 0.9382 0.9392 0.9409 0.9437 0.9437 0.9437 0.9422...
     0.9421 0.9405 0.9360 0.9348 0.9332 0.9283 0.9238 0.9214 0.9199 0.9199 0.9176 0.9108 0.9031 0.9005 0.8971 0.8946 0.8851 0.8702...
     0.84052 0.8272 0.8096 0.7997 0.7849 0.7584 0.7074 0.6527 0.6004 0.4399 0.3540 0.2616 0.2015 0.1671 0.1137 0.05850 0.03400 0 0 0 ...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % External Quantum Efficancy of silicon solar cell

Type_I_Absorption = [0.002 0.00175 0.0015 0.00139 0.00128 0.00117 0.00106 0.00095 0.00084 0.00073 0.00062 0.00052 0.00042 0.00036 0.00032 0.00028...
     0.00026 0.00023 0.00021 0.0002 0.00019 0.00019 0.00018 0.0002 0.00023 0.00027 0.00033 0.0004 0.00047 0.00055 0.00063 0.00084 0.00194 0.00263...
     0.00316 0.00374 0.00403 0.00504 0.00885 0.01274 0.01585 0.02477 0.02629 0.02611 0.02577 0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681...
     0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.2042 0.148 0.1953 0.29512 0.43026 0.65599 1.1591 1.204...
     1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454 5.7883 5.6096 6.8603...
     9.4929 21.630 119.64 110.14 82 65 47 35.5 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601.0...
     337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_IA_Absorption = [0.0022 0.002 0.0018 0.00167 0.00155 0.00142 0.00129 0.00116 0.00104 0.00091 0.00078 0.00068 0.00057 0.00049 0.00044 0.00038 0.00035...
     0.00032 0.0003 0.00028 0.00026 0.00026 0.00025 0.00026 0.00029 0.00032 0.00038 0.00045 0.00052 0.00059 0.00067 0.00089 0.00199 0.00268 0.00322 0.00384 0.00413...
     0.00514 0.00885 0.01274 0.01585 0.02477 0.02629 0.02611 0.02577 0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405...
     0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.2042 0.148 0.1953 0.29512 0.43026 0.65599 1.1591 1.2040 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691...
     31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454 5.7883 5.6096 6.8603 9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592...
     36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601 337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_IB_Absorption = [0.0028 0.0025 0.0022 0.00205 0.0019 0.00175 0.0016 0.00145 0.0013 0.00115 0.001 0.00086 0.00073 0.00063 0.00057 0.00051 0.00047 0.00044 0.00041 0.00038...
     0.00036 0.00035 0.00034 0.00035 0.00038 0.00042 0.00047 0.00052 0.00058 0.00065 0.00072 0.00094 0.00204 0.00273 0.00327 0.00389 0.00418 0.00522 0.00885 0.01274 0.01585...
     0.02477 0.02629 0.02611 0.02577 0.024950 0.022710 0.019640 0.022340 0.027740 0.030910 0.036810 0.048980 0.057940 0.071460 0.092050 0.14405 0.18505 0.23792 0.38759 0.45057...
     0.43851 0.39527 0.20420 0.14800 0.19530 0.29512 0.43026 0.65599 1.1591 1.2040 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691 31.041 25.643 18.833 14.384 11.743 9.9886...
     8.5127 7.8609 6.8674 6.2454 5.7883 5.6096 6.8603 9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401.0...
     6581 3630.8 1402.7 601 337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_II_Absorption = [0.0045 0.0041 0.0037 0.00346 0.00321 0.00297 0.00273 0.00248 0.00224 0.00199 0.00175 0.00154 0.00133 0.00117 0.00106 0.00096 0.0009 0.00084 0.00078 0.00073...
     0.00068 0.00066 0.00063 0.00064 0.00067 0.0007 0.00072 0.00075 0.00079 0.00084 0.00089 0.0011 0.00219 0.0029 0.00348 0.00416 0.00447 0.00552 0.00885 0.01274 0.01585 0.02477...
     0.02629 0.02611 0.02577 0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.20420...
     0.148 0.1953 0.29512 0.43026 0.65599 1.1591 1.2040 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454...
     5.7883 5.6096 6.8603 9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601 337.34...
     163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_III_Absorption = [0.009 0.00775 0.0065 0.00609 0.00568 0.00526 0.00485 0.00444 0.00403 0.00361 0.0032 0.0028 0.0024 0.00213 0.00199 0.00185 0.00175 0.00165 0.00155 0.00145...
     0.00135 0.00127 0.0012 0.00116 0.00115 0.00115 0.00115 0.00116 0.00117 0.00118 0.0012 0.00142 0.00254 0.00327 0.00389 0.00463 0.00499 0.00604 0.00885 0.01274 0.01585 0.02477...
     0.02629 0.02611 0.02577 0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.20420...
     0.148 0.1953 0.29512 0.43026 0.65599 1.1591 1.2040 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454...
     5.7883 5.6096 6.8603 9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601 337.34 163.0...
     122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_1_Absorption = [0.0024 0.0102 0.018 0.01725 0.0165 0.01575 0.015 0.01425 0.0135 0.01275 0.012 0.0104 0.0088 0.00742 0.00626 0.0051 0.0045 0.0039 0.00338 0.00294 0.0025 0.00218...
     0.00186 0.00164 0.00152 0.0014 0.00136 0.00132 0.00128 0.00124 0.0012 0.00144 0.00258 0.00328 0.00386 0.00464 0.00493 0.00594 0.00885 0.01274 0.01585 0.02477 0.02629 0.02611...
     0.02577 0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.2042 0.148 0.1953...
     0.29512 0.43026 0.65599 1.1591 1.2040 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454 5.7883...
     5.6096 6.8603 9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601.0...
     337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_3_Absorption = [0.032 0.028 0.024 0.02313 0.02225 0.02138 0.0205 0.01963 0.01875 0.01788 0.017 0.0146 0.0122 0.01036 0.00908 0.0078 0.00684 0.00588 0.0051 0.0045 0.0039...
     0.0035 0.0031 0.00276 0.00248 0.0022 0.00212 0.00204 0.00198 0.00194 0.0019 0.00206 0.00296 0.00358 0.00412 0.00484 0.00532 0.0065 0.00885 0.01274 0.01585 0.02477 0.02629...
     0.02611 0.02577 0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.2042...
     0.148 0.1953 0.29512 0.43026 0.65599 1.1591 1.2040 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674...
     6.2454 5.7883 5.6096 6.8603 9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601.0...
     337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_5_Absorption = [0.04 0.0375 0.035 0.0335 0.032 0.0305 0.029 0.0275 0.026 0.0245 0.023 0.0202 0.0174 0.015 0.013 0.011 0.00972 0.00844 0.00736 0.00648 0.0056 0.00508 0.00456...
     0.00416 0.00388 0.0036 0.0034 0.0032 0.00308 0.00304 0.003 0.00324 0.0038 0.00432 0.00492 0.00566 0.00619 0.0074 0.00885 0.01274 0.01585 0.02477 0.02629 0.02611 0.02577...
     0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.2042 0.148 0.1953 0.29512...
     0.43026 0.65599 1.1591 1.2040 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.5830 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454 5.7883 5.6096...
     6.8603 9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601.0...
     337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_7_Absorption = [0.048 0.044 0.04 0.03875 0.0375 0.03625 0.035 0.03375 0.0325 0.03125 0.03 0.0264 0.0228 0.02 0.018 0.016 0.0144 0.0128 0.01138 0.01014 0.0089 0.00818 0.00746...
     0.00684 0.00632 0.0058 0.00544 0.00508 0.00484 0.00472 0.0046 0.00460 0.00474 0.00504 0.00558 0.00666 0.00738 0.00864 0.00885 0.01274 0.01585 0.02477 0.02629 0.02611 0.02577...
     0.02495 0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.2042 0.148 0.1953 0.29512...
     0.43026 0.65599 1.1591 1.204 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.583 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454 5.7883 5.6096 6.8603...
     9.4929 21.630 119.64 110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601.0...
     337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

Type_9_Absorption = [0.05 0.046 0.042 0.04163 0.04125 0.04088 0.0405 0.04013 0.03975 0.03938 0.039 0.0354 0.0318 0.0288 0.0264 0.024 0.022 0.02 0.0184 0.0172 0.016 0.01452 0.01304...
     0.01182 0.01086 0.0099 0.00906 0.00822 0.0075 0.0069 0.0063 0.0059 0.00594 0.0062 0.00672 0.00798 0.00875 0.01028 0.00885 0.01274 0.01585 0.02477 0.02629 0.02611 0.02577 0.02495...
     0.02271 0.01964 0.02234 0.02774 0.03091 0.03681 0.04898 0.05794 0.07146 0.09205 0.14405 0.18505 0.23792 0.38759 0.45057 0.43851 0.39527 0.2042 0.148 0.1953 0.29512 0.43026 0.65599...
     1.1591 1.204 1.2566 1.1552 1.0777 1.1885 1.8188 3.7699 9.583 31.691 31.041 25.643 18.833 14.384 11.743 9.9886 8.5127 7.8609 6.8674 6.2454 5.7883 5.6096 6.8603 9.4929 21.630 119.64...
     110.14 82 65 47 35.500 27.797 21.675 19.318 22.592 36.485 69.169 95.504 122 675.83 5432.4 11325 11401 6581 3630.8 1402.7 601.0...
     337.34 163 122.20 112.44 122.47 144.55]; %(cm^(-1))

% Physical constants
e = 1.602*10^-19; %Charge of electron [Coulombs or J V^(-1)]
m = 1.5; %Operating voltage correction factor; ASSUMPTION
k = 1.381*10^-23; %Boltzmanns's constant [m^2 kg s^-2 K^-1 or J K^(-1)]
c = 299792458; %speed of light [m s^(-1)]
h = 6.62607004*10^(-34); %Plancks constant [m^(2) kg s^(-1)]

% I_o = 2.36*10^-8; %Saturation current for silicon [A/cm^2]; ASSUMPTION
% I_o = 1*10^-12; %Saturation current for silicon [A/cm^2]; ASSUMPTION GREG ORIGINAL

if Water_Quality == 1
    Absorption = Type_I_Absorption;
elseif Water_Quality == 2
    Absorption = Type_IA_Absorption;
elseif Water_Quality == 3
    Absorption = Type_IB_Absorption;
elseif Water_Quality == 4
    Absorption = Type_II_Absorption;
elseif Water_Quality == 5
    Absorption = Type_III_Absorption;
elseif Water_Quality == 6
    Absorption = Type_1_Absorption;
elseif Water_Quality == 7
    Absorption = Type_3_Absorption;
elseif Water_Quality == 8
    Absorption = Type_5_Absorption;
elseif Water_Quality == 9
    Absorption = Type_7_Absorption;
elseif Water_Quality == 10
    Absorption = Type_9_Absorption;
end

%% Calculations
% IXYS Characteristics
SR = (e.*EQE.*Lambda)./(h.*c.*10^6); %Spectral response [A W^(-1)]; 10^6 converts c to [microm/s]
Eg = (Eg_0 - ((alpha.*T.^2)./(T+beta))).*1.60218.*10^(-19); %band gap energy [J]
Jo = C.*T.^(3).*exp(-Eg./(m.*k.*T)); %saturation current density [A m^(-2)]
Io = Jo.*Solar_Cell_Area;

% Power
Solar_Spectrum_Incident_To_Panel = (Scaled_Direct_Irradiance_Spectrum.*repmat(cosd(Incidence_Angle'),122,1)+Scaled_Diffuse_Irradiance_Spectrum.*repmat(cosd(Panel_Tilt_Bound'./2),122,1).^2).*(repmat(Depth',122,1) == 0) + ... %Solar irradiance normal to PV cell when surfaced
    (Scaled_Transmitted_GHI_Spectrum.*exp(-repmat(Absorption',1,length(Depth)).*repmat(Depth',122,1).*100).*repmat(cosd(Panel_Tilt_Bound'),122,1)).*(repmat(Depth',122,1) ~= 0); %Solar irradiance normal to PV cell when submerged [W m^(-2)]
JL_Spectrum = repmat(SR',1,length(Depth)).*Solar_Spectrum_Incident_To_Panel; %Light current denisty at each wavelength [A m^(-2)]
JL = trapz(Lambda,JL_Spectrum)'; %Light current denisty [A m^(-2)]
IL = JL.*Solar_Cell_Area; %cell light current [A]
V_oc = ((m.*k.*T)./e).*log((IL./Io)+1); %cell open circuit voltage [V]

%set physical limit of IXYS Si cell
IL(IL > 0.025) = 0.025; 
% V_oc(V_oc > 6.3) = 6.3; %I don't think we need to limit voc because IL v
% Io relationship shouldn't let it reach too high

%account for array effects on voltage and current
switch array_type
    case 'series'
        V_oc = V_oc.*array_n; %when in series voltages add 
    case 'parallel'
        IL = IL.*array_n; %current produced by one cell is I_total/array_n so I_total = Isc*array_n
end

Rch = V_oc./IL;   %characteristic resistence [ohm]
rs = Rs./Rch;   %normalized series resistence [-]
v_oc = (V_oc.*e)./(m.*k.*T);    %normalized open cicuit voltage [-]
I_n = 0:0.01:1;
fi = repmat(I_n,length(Depth),1);  %normalized current, will range from 0 to 1 for all curves (regardless of light intensity)
FFs = max(fi.*((1-fi.*repmat(rs,1,length(I_n)))+(1./repmat(v_oc,1,length(I_n))).*log(1-fi.*(1-exp(repmat(v_oc,1,length(I_n)).*repmat(rs,1,length(I_n))-1)))),[],2);    %FF for ignoring the effects of Rsh, but accounting for Rs; not restricted by v_oc > 10
%when there is no cell illumination
%if V_oc & I_sc = 0    %if V_oc = 0           %if I_sc = 0
% R_ch = nan           % R_ch = Inf           % R_ch = 0
% rs = nan             % rs = 0               % rs = Inf
% v_oc = 0             % v_oc = number        % v_oc = 0
% FFs = nan            % FFs = number         % FFs = nan
Temp_FF = (V_oc == 0) + (IL ==0);
FFs(Temp_FF >= 1) = 0;    %this will remove imaginary values caused by 1/0 in FFs due to 0 voltage or current output
Solar_Cell_Power_Output = FFs.*IL.*V_oc;  %Power Output [W]

Total_Solar_Irradiance_Incident_To_Cell = trapz(Lambda,Solar_Spectrum_Incident_To_Panel)';
end