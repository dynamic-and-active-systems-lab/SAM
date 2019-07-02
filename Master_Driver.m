% close all
clear all
clc

%% Model parameters/User Interaction
%{
The code will run through a series of option for the user's upload. It can
read in a single file containing the deployment data. The current supported 
data files are Excel, CSV, and MAT files.

If formated as an .xlsx and or .csv file, the input data should be compiled
where each column represents:
Time | Depth | Latitude | Longitude | Temperature | Ax | Ay | Az |
If formatted as a mat file, each variable named as follows:
time | depth |   lat    |   lon     |    temp     | ax | ay | az |
If a user does not have a certain variable the respective column should be
left blank, but the order/spacing should be maintained.
%}
%File selection
%must be in current directory (GUI should have browser functionality)
in_filepath = 'C:\Users\amysw\OneDrive - Northern Arizona University\Marine Energy Harvesting\';
% in_filepath = 'C:\Users\ams883\OneDrive - Northern Arizona University\Marine Energy Harvesting\';%'.';
file = strcat(in_filepath,'2018007_17A0036_Kilo.mat');%input('Please input your filename.ext: ','s');
out_filepath = strcat(in_filepath,'Output');%'./Output';
out_filetpe = 'mat';
% C = strsplit(file,"."); 
base = '2018007_17A0036_Kilo';%string(C(1)); 
ext = 'mat';%string(C(2)); %split up the input filename into the base and extension

% Water quality
struct2cell_str = 'ocean';
switch struct2cell_str % Ranges 1-10. 1-5 open ocean, 6-10 coastal
    case 'ocean'
        Water_Quality = [1 5];
    case 'coast'
        Water_Quality = [6 10];
    case 'custom'
        Water_Quality = input('Choose your water quality/qualities. Ranges 1-10; 1-5 open ocean, 6-10 coastal\n');
end
%Time zone
Time_Zone = -7;%input('What time zone was the data collected in? (UTC + X): '); %-7;  %UTC + X 

%CERES data folder location
% irr_file_path = 'C:\Users\laa22\Desktop\CERES Irradiance Data'; %laptop path
irr_filepath = strcat(in_filepath,'CERES Data 2018');%'C:\Users\laa22\OneDrive\Desktop\CERES Irradiance Data'; %home computer path
% irr_file_path = 'C:\Users\laa222\Desktop\CERES Irradiance Data'; %school computer path
% irr_file_path = '\\EGRSHARES\Homes\NAU\laa222\Desktop\CERES Irradiance Data'; %remote desktop path

% IXYS Characteristic values (suggest these values as default)
array_type = 'series'; %can be series or parallel
array_n = 10; %number of elements in the array
Solar_Cell_Area = 0.01.*0.006; %units in m^2 (10mm x 6mm)
Rsh = 50000; %shunt resistence of cell; empirically derived in lab experiments [ohm]
Rs = 4.82;   %series resistence of cell empirically derived from in lab experiments [ohm]

%%%%%%%%%% ADVANCED %%%%%%%%%%
% Saturation Current constants (
%for Si
C = 500; %current density constant [A m^(-2) K^(-3)] 
Eg_0 = 1.1557; %[eV]
alpha = 7.021*10^(-4); %[eV K^(-1)]
beta = 1108; %[K]
%% Data Formatting
%{
NOTE: This function previously supported reading in multiple files,
however, this was removed for ease of maintenance. It now only reads in 
single files, as described above. If the previous version of this function 
is desired, itcan be found in the '2018-06 - Open Sourcing Robustness'
folder in the 'Model Edits'
%}
[Year, Month, Day] = Data_Format_Interpolation(file, base, ext, out_filetpe); %pass the read in function the relevant information

%% Outline Logic for Month Run
% outputNames = {'Zenith', 'Total_Clear_Sky_Direct', 'Total_Clear_Sky_GHI', 'CERES_Direct_Down', 'CERES_GHI', 'Total_Scaled_GHI', 'Total_Scaled_Transmitted_GHI', 'Total_Solar_Irradiance_Incident_To_Cell', 'Solar_Cell_Power_Output', 'V_oc', 'I_sc'}; %would be defined by the user through the GUI interface
outputNames = {'Time', 'Depth', 'Latitude', 'Longitude', 'T', 'Panel_Tilt', 'Panel_Azimuth', 'Incidence_Angle', 'Solar_Cell_Power_Output', 'Total_Scaled_GHI', 'Scaled_GHI_Spectrum', 'Total_Solar_Irradiance_Incident_To_Cell', 'V_oc', 'I_sc'};

wb_progress = 0; %varible to indicate percent processing complete
wb_unit = 1/(length(Day)*length(Water_Quality)); %unit of waitbar increment
wb = waitbar(wb_progress,'Calculating');
filename = string(ls(sprintf('%s/%s_*.mat', in_filepath, base))); 
for iWQ = 1:length(Water_Quality)
    for jfile = 1:length(filename)
        wb_str = sprintf('Water Quality %d/%d | Month %d/%d \n',iWQ,length(Water_Quality),jfile,length(filename)); %string concerning what month of deployment are being conducted as MAIN.m does not keep track of this information
        [out, wb_progress] = MAIN(wb, wb_progress, wb_str, wb_unit, filename(jfile), Water_Quality(iWQ), Time_Zone, outputNames, in_filepath, irr_filepath, array_type, array_n, Solar_Cell_Area, Rsh, Rs, C, Eg_0, alpha, beta); %pass to main model
        eval(sprintf('out_%d_WQ%d=out;',jfile,Water_Quality(iWQ))) %save main model output as unique structure name for later identification       
        clear out
    end

%% Concatenate Output Structure Variables
%{
the cellfun function is able to concatenate structures taht share the same 
fieldnames, the function should look like:
S = cell2struct(cellfun(@horzcat,struct2cell(Struct1),struct2cell(Struct2),...,struct2cell(SN),'UniformOutput',false),fieldnames(Struct1),1)
where S is the new concatenated structure, however, since there will be
a variable length (deployments are variable months long), we use a for
loop to cycle through the outputs of each month of the deployment and
create the struct2cell(Struct1),...,struct2cell(StructN) portion of the
statement.
%}

waitbar(0,wb,'Saving Output') %reset the waitbar

waitbar((iWQ-1)/length(iWQ),wb,'Saving Output') %iterate the waitbar

struct2cell_str = ''; %reset struct2cell string for cellfun function
structNames = who(sprintf('out_*WQ%d',Water_Quality(iWQ))) %find month's output structure
for istruct = 1:length(structNames) %for each month output structure
    struct2cell_str = [struct2cell_str,sprintf('struct2cell(%s),',structNames{istruct})]; %print struct2cell(struct1),...,struct2cell(struct2),... etc.
end
%assemble full cellfun function, as described above
cellfun_str = sprintf('out_WQ%d=cell2struct([%s],fieldnames(%s),1);',Water_Quality(iWQ),struct2cell_str,structNames{1}) %concatenate structures
eval(cellfun_str); %evaluate the function

%calculate Energy of deployment; since it is calculated using a
%cumulative summation of power it must be calculated after
%concatenation
eval(sprintf('temp = out_WQ%d.Time;',Water_Quality(iWQ)))
energy_str = sprintf('out_WQ%d.Solar_Cell_Energy_Output = reshape(cumtrapz(datenum(temp),reshape(out_WQ%d.Solar_Cell_Power_Output,1,[]).*86400),86400,[]);',Water_Quality(iWQ),Water_Quality(iWQ))
eval(energy_str)

%save
if isdir(out_filepath) == 0; mkdir(out_filepath); end % if a directory for the output doesn't exist, create one
save(sprintf('%s/%s_WQ%dOutput.mat',out_filepath, base, Water_Quality(iWQ)),'-struct',sprintf('out_WQ%d',Water_Quality(iWQ)))
clearvars -except in_filepath out_filepath irr_filepath base ext Year Month Day Water_Quality Time_Zone outputNames wb_progress wb_unit wb filename iWQ jfile wb_str array_type array_n Solar_Cell_Area Rsh Rs C Eg_0 alpha beta
end
close(wb)

% toc