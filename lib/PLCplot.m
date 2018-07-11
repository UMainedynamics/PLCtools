% These functions allow for the use of the Data Cursor to pull information
% from the plc object and plot it alongside its normalized coordinates.

function plothandle = PLCplot(fighandle,plcobj,plcData)
% Plots graph and sets up a custom data tip update function
% for a plot with one component
plothandle = pcolor(plcobj.xgrid,plcobj.ygrid,plcData);
plothandle.EdgeColor = 'none';
dcm_obj = datacursormode(fighandle);
set(dcm_obj,'UpdateFcn',{@PlotDataCursorText,plcobj,plcData})
end

function txt = PlotDataCursorText(~,event_obj,plcobj,plcData)
% Customizes text of data tips
pos = get(event_obj,'Position');
txt = {...
    ['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))],...
    ['Value: ', sprintf('%.2e',plcData( ...
        round((pos(2)/max(plcobj.total_size))*plcobj.msResolution), ...
        round((pos(1)/max(plcobj.total_size))*plcobj.msResolution)))]};
end