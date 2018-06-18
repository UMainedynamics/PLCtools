cd('C:\Users\NWRichmond\Documents\GitHub\PLCtools\DEMO_DATA\Fold3cgsm30_2m_9d_11h_59m')
load('Fold3cgsm30.mat')
mst =[ms.MicroCoordinates, ms.MicroStress{10, 1}(:,end)];
mst_table = array2table(mst,...
    'VariableNames',{'x','y','stress'});
F = scatteredInterpolant(mst_table.x,mst_table.y,mst_table.stress);
[xq,yq] = meshgrid(linspace(0,1,100),linspace(0,0.75,100));
cq = F(xq,yq);
h = pcolor(xq,yq,cq);
h.EdgeColor = 'none';

figure
phase_to_exclude = 2;
[row,~] = find(ms.GrainPhases == phase_to_exclude);
scatter(ms.Grains{1, row}(:,1),ms.Grains{1, row}(:,2))

% hold on
% scatter(ms.Grains{1, 1}(:,1),ms.Grains{1, 1}(:,2))
% scatter(ms.Grains{1, 2}(:,1),ms.Grains{1, 2}(:,2))

% scatter(ms.GrainPolylines{row, 1}(:,1),ms.GrainPolylines{row, 1}(:,2))
% hold on
% scatter(ms.GrainPolylines{row, 2}(:,1),ms.GrainPolylines{row, 2}(:,2))
% scatter(ms.GrainPolylines{row, 3}(:,1),ms.GrainPolylines{row, 3}(:,2))
% scatter(ms.GrainPolylines{row, 4}(:,1),ms.GrainPolylines{row, 4}(:,2))