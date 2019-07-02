function [Output_Name,Units] = Output_var2name(Output_Variable)
%Find full name of an output variable
    % Input: cell array of strings

for i = 1:length(Output_Variable)
    switch Output_Variable{i}
        case 'Time'
            Output_Name{i} = 'Time';
            Units{i} = Output_Name{i};
        case 'Depth'
            Output_Name{i} = 'Depth';
            Units{i} = strcat(Output_Name{i},' (m)');
        case 'Latitude'
            Output_Name{i} = 'Latitude';
            Units{i} = Output_Name{i};
        case 'Longitude'
            Output_Name{i} = 'Longitude';
            Units{i} = Output_Name{i};
        case 'T'
            Output_Name{i} = 'Temperature';
            Units{i} = strcat(Output_Name{i},' ($^\circ$C)');
        case 'Panel_Tilt'
            Output_Name{i} = 'Panel Tilt';
            Units{i} = strcat(Output_Name{i},' ($^\circ$)');
        case 'Panel_Azimuth'
            Output_Name{i} = 'Panel Azimuth';
            Units{i} = strcat(Output_Name{i},' ($^\circ$)');
        case 'Zenith'
            Output_Name{i} = 'Panel Zenith';
            Units{i} = strcat(Output_Name{i},' ($^\circ$)');
        case 'Incidence_Angle'
            Output_Name{i} = 'Angle of Incidence';
            Units{i} = strcat(Output_Name{i},' ($^\circ$)');
        case 'Solar_Cell_Power_Output'
            Output_Name{i} = 'Maximum Power';
            Units{i} = strcat(Output_Name{i},' (W)');
        case 'Solar_Cell_Energy_Output'
            Output_Name{i} = 'Maximum Energy';
            Units{i} = strcat(Output_Name{i},' (J)');
        case 'Total_Scaled_GHI'
            Output_Name{i} = 'Global Horizantal Irradiance (Total)';
            Units{i} = strcat(Output_Name{i},' (W/m$^2$)');
        case 'Scaled_GHI_Spectrum'
            Output_Name{i} = 'Global Horizantal Irradiance (Spectrum)';
            Units{i} = strcat(Output_Name{i},' (W/m$^2$)');
        case 'Total_Solar_Irradiance_Incident_To_Cell'
            Output_Name{i} = 'Total Irradiance Incident to Cell';
            Units{i} = strcat(Output_Name{i},' (W/m$^2$)');
        case 'V_oc'
            Output_Name{i} = 'Open Current Voltage';
            Units{i} = strcat(Output_Name{i},' (V)');
        case 'I_sc'
            Output_Name{i} = 'Light Current';
            Units{i} = strcat(Output_Name{i},' (A)');
    end
end
end

