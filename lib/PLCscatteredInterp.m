user = 'NWRichmond'; % 'Peter Koons' for work pc or 'NWRichmond' for home pc
plc_case = 'Fold3cgsm30_2m_9d_11h_59m';
cd(strcat('C:\Users\', ...
    user,'\Documents\GitHub\PLCtools\', ...
    'DEMO_DATA\',plc_case));
load('Fold3cgsm30.mat')
resolution = 100;
phase_to_exclude = 2;
[row,~] = find(ms.GrainPhases == phase_to_exclude);
[row_elem,~] = find(ms.ElementPhases == phase_to_exclude);

coords = ms.MicroCoordinates/max(ms.MicroCoordinates(:));
mst =[coords, ms.MicroStress{10, 1}(:,end)];
mst_table = array2table(mst, ...
    'VariableNames',{'x','y','stress'});
mst_mask = mst(row_elem,:);
mst_mask_table = array2table(mst_mask,...
    'VariableNames',{'x','y','stress'});

F = scatteredInterpolant(mst_table.x,mst_table.y,mst_table.stress);
max_x = max(mst_table.x);
max_y = max(mst_table.y);
[xq,yq] = meshgrid(linspace(0,max_x,max_x*resolution), ...
    linspace(0,max_y,max_x*resolution));
cq = F(xq,yq);

xGrain = [ms.Grains{1, row}(:,1),ms.Grains{1, row}(:,2)];
xGrainNorm = xGrain/max(xGrain(:));
xv = xGrainNorm(:,1);
yv = xGrainNorm(:,2);
[in,on] = inpolygon(xq,yq,xv,yv);
cq(in) = NaN;
cq(on) = NaN;

h = pcolor(xq,yq,gradient(cq));
h.EdgeColor = 'none';

% hold on
% px = polyshape(xGrainNorm,'Simplify',true);
% plot(px)


% scatter(mst_mask(:,1),mst_mask(:,2))
% scatter(xGrainNorm(:,1),xGrainNorm(:,2))

% hold on
% scatter(ms.Grains{1, 1}(:,1),ms.Grains{1, 1}(:,2))
% scatter(ms.Grains{1, 2}(:,1),ms.Grains{1, 2}(:,2))

% scatter(ms.GrainPolylines{row, 1}(:,1),ms.GrainPolylines{row, 1}(:,2))
% hold on
% scatter(ms.GrainPolylines{row, 2}(:,1),ms.GrainPolylines{row, 2}(:,2))
% scatter(ms.GrainPolylines{row, 3}(:,1),ms.GrainPolylines{row, 3}(:,2))
% scatter(ms.GrainPolylines{row, 4}(:,1),ms.GrainPolylines{row, 4}(:,2))