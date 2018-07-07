cutline = auto_xc_linelength;
if size(cutline) == 1
    obj.CutLineMaxGradient(cutline);
elseif size(cutline) == 2
    verts = cutline;
    obj.vx1 = obj.Denormalize(verts(1,1));
    obj.vy1 = obj.Denormalize(verts(1,2));
    obj.vx2 = obj.Denormalize(verts(2,1));
    obj.vy2 = obj.Denormalize(verts(2,2));
end
    obj.xc_x = [obj.vx1,obj.vx2];
    obj.xc_y = [obj.vy1,obj.vy2];
    obj.cut_line_v1 = [obj.vx1,obj.vy1]; % cartesian coordinates for
        % the 1st vertex of the cut line / line of cross section
    obj.cut_line_v2 = [obj.vx2,obj.vy2]; % cartesian coordinates for
        % the 2nd vertex of the cut line / line of cross section
    obj.vertex_distance = round(pdist([obj.vx1,obj.vy1;...
        obj.vx2,obj.vy2],'euclidean'));
    obj.cutlinelength = obj.vertex_distance;
    xc_line = smooth(improfile(obj.gridded_data,...
        [obj.vx1,obj.vx2], [obj.vy1,obj.vy2], obj.vertex_distance+1,...
        'bilinear'),'moving');  % find the values along the line
    obj.cutline = xc_line;