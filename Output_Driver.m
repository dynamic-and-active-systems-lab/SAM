function Output_Driver(in_filepath, base, out_filepath, irr_filepath, outputNames,Day, Water_Quality, Time_Zone, array_type, array_n, Solar_Cell_Area, Rsh, Rs, C, Eg_0, alpha, beta,Cust_Ext_Eff)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
wb_progress = 0; %varible to indicate percent processing complete
wb_unit = 1/(length(Day)*length(Water_Quality)); %unit of waitbar increment
wb = waitbar(wb_progress,'Calculating');
filename = string(ls(sprintf('%s/%s_*.mat', in_filepath, base))); 
for iWQ = 1:length(Water_Quality)
    % Check if output files already exist (function takes a long time to run, so the user might have run it partially)
    out_filename = sprintf('%s/%s_WQ%dOutput.mat',out_filepath, base, Water_Quality(iWQ));
    if isfile(out_filename)==0 % file doesn't exist, continue with calcs
        for jfile = 1:length(filename)
            wb_str = sprintf('Water Quality %d/%d | Month %d/%d \n',iWQ,length(Water_Quality),jfile,length(filename)); %string concerning what month of deployment are being conducted as MAIN.m does not keep track of this information
            [out, wb_progress] = MAIN(wb, wb_progress, wb_str, wb_unit, filename(jfile), Water_Quality(iWQ), Time_Zone, outputNames, in_filepath, irr_filepath, array_type, array_n, Solar_Cell_Area, Rsh, Rs, C, Eg_0, alpha, beta,Cust_Ext_Eff); %pass to main model
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

%         struct2cell_str = ''; %reset struct2cell string for cellfun function
        struct2cell_cell = {};
        structNames = who(sprintf('out_*WQ%d',Water_Quality(iWQ))); %find month's output structure
        for istruct = 1:length(structNames) %for each month output structure
%             struct2cell_str = [struct2cell_str,sprintf('struct2cell(%s),',structNames{istruct})]; %print struct2cell(struct1),...,struct2cell(struct2),... etc.
            struct2cell_cell = [struct2cell_cell,struct2cell(eval(structNames{istruct}))];
        end
        %assemble full cellfun function, as described above
        out_WQ = cell2struct(struct2cell_cell,fieldnames(eval(structNames{1})),1);
%         assignin('base','out_WQ',out_WQ)
%         whos out_WQ
%         cellfun_str = [sprintf('out_WQ%d',Water_Quality(iWQ)),sprintf('=cell2struct([%s],fieldnames(%s),1);',struct2cell_str,structNames{1})] %concatenate structures
%         eval(cellfun_str); %evaluate the function
    
        %calculate Energy of deployment; since it is calculated using a
        %cumulative summation of power it must be calculated after concatenation
        time_temp = out_WQ.Time; %whos time_temp
        power_temp = out_WQ.Solar_Cell_Power_Output; %whos power_temp
        out_WQ(1).Solar_Cell_Energy_Output = reshape(cumtrapz(datenum(time_temp),reshape(power_temp,1,[]).*86400),86400,[]);
        % Old method (doesn't work for new output variables)
%         energy_str = [sprintf('out_WQ%d.Solar_Cell_Energy_Output',Water_Quality(iWQ)),sprintf('= reshape(cumtrapz(datenum(out_WQ%d.Time),reshape(out_WQ%d.Solar_Cell_Power_Output,1,[]).*86400),86400,[]);',Water_Quality(iWQ),Water_Quality(iWQ))];
%         eval(energy_str)
        

        %save
        if isdir(out_filepath) == 0; mkdir(out_filepath); end % if a directory for the output doesn't exist, create one
            save(out_filename,'out_WQ')
            % issues saving the output as individual variables instead of the structure.  This is addressed when the data is loaded.
%             save(out_filename,'-struct','out_WQ')
%             save(out_filename,'-struct',sprintf('out_WQ%d',Water_Quality(iWQ)))
            clearvars -except in_filepath out_filepath irr_filepath base ext Year Month Day Water_Quality Time_Zone outputNames wb_progress wb_unit wb filename iWQ jfile wb_str array_type array_n Solar_Cell_Area Rsh Rs C Eg_0 alpha beta
        end
    end
close(wb)
end

