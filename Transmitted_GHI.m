function [Total_Scaled_GHI, Scaled_Transmitted_GHI_Spectrum, Total_Scaled_Transmitted_GHI] = Transmitted_GHI( Scaled_Direct_Irradiance_Spectrum,Scaled_Diffuse_Irradiance_Spectrum,Scaled_GHI_Spectrum,Zenith)
% This function calculates the Total GHI and transmitted GHI spectrum

%% Constants

Lambda = [0.300 0.305 0.310 0.315 0.320 0.325 0.330 0.335 0.340 0.345 0.350 0.360 0.370 0.380 0.390 0.400 0.410 0.420 0.430 0.440...
    0.450 0.460 0.470 0.480 0.490 0.500 0.510 0.520 0.530 0.540 0.550 0.570 0.593 0.610 0.630 0.656 0.668 0.690 0.710 0.718...
    0.724 0.740 0.753 0.758 0.763 0.768 0.780 0.800 0.816 0.824 0.832 0.840 0.860 0.880 0.905 0.915 0.925 0.930 0.937 0.948...
    0.965 0.980 0.994 1.040 1.070 1.100 1.120 1.130 1.145 1.161 1.170 1.200 1.240 1.270 1.290 1.320 1.350 1.395 1.443 1.463...
    1.477 1.497 1.520 1.539 1.558 1.578 1.592 1.610 1.630 1.646 1.678 1.740 1.800 1.860 1.920 1.960 1.985 2.005 2.035 2.065...
    2.100 2.148 2.198 2.270 2.360 2.450 2.500 2.600 2.700 2.800 2.900 3.000 3.100 3.200 3.300 3.400 3.500 3.600 3.700 3.800...
    3.900 4.000]; %Wavelengths

N1 = 1.00; %Refractive index for air

N2 = [1.37 1.37 1.37 1.37 1.37 1.36 1.36 1.36 1.36 1.36 1.36 1.36 1.35 1.35 1.35 1.35 1.35 1.35 1.35 1.34 1.34 1.34 1.34 1.34 1.34...
    1.34 1.34 1.34 1.34 1.34 1.34 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33 1.33...
    1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32 1.32...
    1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.30 1.30 1.30 1.30 1.30 1.30 1.30 1.29 1.29 1.29 1.29...
    1.28 1.27 1.26 1.25 1.23 1.18 1.13 1.20 1.35 1.45 1.46 1.43 1.40 1.38 1.37 1.36 1.35 1.34 1.33]; N2 = repmat(N2',1,length(Zenith)); %Refractive Index of Water

%% Percent transmitted into water

Rp = abs((N1.*sqrt(1-(N1./N2.*repmat(sind(Zenith'),length(Lambda),1)).^2)-N2.*repmat(cosd(Zenith'),length(Lambda),1))./(N1.*sqrt(1-(N1./N2.*repmat(sind(Zenith'),length(Lambda),1)).^2)+N2.*repmat(cosd(Zenith'),length(Lambda),1))).^2; %Reflectivity of P Polarized Light
Rs = abs((N1.*repmat(cosd(Zenith'),length(Lambda),1)-N2.*sqrt(1-(N1./N2.*repmat(sind(Zenith'),length(Lambda),1)).^2))./(N1.*repmat(cosd(Zenith'),length(Lambda),1)+N2.*sqrt(1-(N1./N2.*repmat(sind(Zenith'),length(Lambda),1)).^2))).^2;  %Reflectivity of S Polarized Light
Percent_Transmitted = (1-(Rp+Rs)./2); %Percent of irradiance transmitted into water

%% Scaled transmitted GHI

Reflected_Zenith = asind((N1.*repmat(sind(Zenith'),length(Lambda),1))./N2);
Scaled_Transmitted_GHI_Spectrum = Scaled_Direct_Irradiance_Spectrum.*Percent_Transmitted.*cosd(Reflected_Zenith) + Scaled_Diffuse_Irradiance_Spectrum.*.934;
Total_Scaled_Transmitted_GHI = trapz(Lambda,Scaled_Transmitted_GHI_Spectrum(:,:))';

%% Scaled GHI

Total_Scaled_GHI = trapz(Lambda,Scaled_GHI_Spectrum(:,:))';

end