function subhandle = PLCsubplot(fighandle,m,n,p,plcobj,params)
% Plots graph and sets up a custom data tip update function
% for a subplot (either horizontal or vertical) with two components
subhandle(p) = subplot(m,n,p);
parameter = params{1,p};
imagesc(parameter);
dcm_obj = datacursormode(fighandle);
set(dcm_obj,'UpdateFcn',{@SubPlotDataCursorText,plcobj,subhandle(p),p,params})
end

function txt = SubPlotDataCursorText(~,event_obj,plcobj,subhandle,p,params)
% Customizes text of data tips
obj_size_x = plcobj.total_size(1);
obj_size_y = plcobj.total_size(2);
pos = get(event_obj,'Position');
cTarget = get(event_obj.Target,'parent');
if cTarget == subhandle
    parameter = params{1,p};
else
    parameter = params{1,p-1};
end
txt = {...
    ['X: ',num2str(round(pos(1))/obj_size_x)],...
    ['Y: ',num2str(1-(round(pos(2))/obj_size_y))],...
    ['Value: ',num2str(parameter(pos(1),pos(2)))]};
end