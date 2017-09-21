% AUTHOR: Nick Richmond
% Created September 2017
%
% DESCRIPTION: Creates files for batch analysis for the Power-law Creep
% Toolbox (PLC)

%% SETUP

picspath = fullfile(userpath,'PLCToolboxV06_64bit',...
    'Custom and Batch Analyses','Sawtooth','Frequency'); 
    % use pwd for current directory
project_basename = 'Sawtooth_frequency'; % basename for pictures
load_type = 'stress'; % either 'stress' or 'strain rate'
phase_names = {'weak','strong'};
Ad = {1e-21,1e-19};
n = {3,3};
Q = {0,0};
mesh_density = 0.025; % between 0.005 and 0.2
temperature(1).value = 773;
temperature(1).tensor(1).value = [0;0;0;0;0;5e9]; %[11,22,33,32,31,12)
temperature(1).tensor(2).value = [0;0;0;0;0;6e9];
%%% EXAMPLE: ADDING A SECOND TEMPERATURE
% temperature(2).value = 823;
% temperature(2).tensor(1).value = [0;0;0;0;0;5.5e9]; %[11,22,33,32,31,12)
% temperature(2).tensor(2).value = [0;0;0;0;0;6.5e9];

%% WRITE THE ANALYSIS FILES

if strcmp(picspath,pwd)~=1
    cd(picspath);
end
if strcmp(load_type,'stress')==1
    load_type_value = 0;
elseif strcmp(load_type,'strain rate')==1
    load_type_value = 1;
end

pics = dir(strcat(project_basename,'*.png'));
batch_size = length(pics);
fname = cell(batch_size,1);
for ia = 1:batch_size
    [~,picname,~] = fileparts(pics(ia).name);
    fname{ia} = strcat(picname,'_runfile','.txt');
    fid = fopen(fname{ia},'w');
    fprintf(fid,'%s\r\n','#name of image file',pics(ia).name,newline);
    fprintf(fid,'%s\r\n','#Load type',num2str(load_type_value),newline);
    fprintf(fid,'%s\r\n','#Number of phases',num2str(length(phase_names)),newline);
    for ib = 1:length(phase_names)
        phase_name_spec = '#phase %d name, Ad, n, and Q values \r\n';
        fprintf(fid, phase_name_spec,ib);    
        fprintf(fid,'%s\r\n',num2str(phase_names{ib}),num2str(Ad{ib}),...
            num2str(n{ib}),num2str(Q{ib}),newline);
    end
    fprintf(fid,'%s\r\n','#Number of Temperature Intervals',...
        num2str(length(temperature)),newline);
    for ic = 1:length(temperature)
        fprintf(fid,'%s\r\n',['#Temperature ',num2str(ic)],...
            num2str(temperature(ic).value),newline);
        fprintf(fid,'%s\r\n',['#Number of ',load_type,...
            ' events for temperature #',num2str(ic)],...
            num2str(length(phase_names)),newline);
        for id = 1:length(temperature(ic).tensor)
            fprintf(fid,'%s\r\n',['#',load_type,...
                ' components for temperature #',num2str(ic)]);
                fprintf(fid,'%#.2e\r\n',...
                    temperature(ic).tensor(id).value(1),...
                    temperature(ic).tensor(id).value(2),...
                    temperature(ic).tensor(id).value(3),...
                    temperature(ic).tensor(id).value(4),...
                    temperature(ic).tensor(id).value(5),...
                    temperature(ic).tensor(id).value(6),[ ],[ ]);
        end
    end
    fprintf(fid,'%s\r\n','#mesh density value (between 0.005 and 0.2)');
    fprintf(fid,'%g\r\n',mesh_density);
    fclose(fid);
end

%% GENERATE THE BATCH CONTROL FILE

ctrl = fopen([project_basename,'_Batch','.txt'],'w');
for ie = 1:batch_size
    fprintf(ctrl,'%s\r\n',fname{ie});
end
fclose(ctrl);
