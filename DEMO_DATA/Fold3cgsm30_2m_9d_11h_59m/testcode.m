[maxgradrow,maxgradcol] = find(plc.parameter_gradient == plc.max_gradient);
plc.max_gradient_row = maxgradrow;
plc.max_gradient_col = maxgradcol;
plc.vxcenter = maxgradcol;
plc.vycenter = maxgradrow;


%%% GRADIENT LINE ROTATION
    alpha = pi/144; % angle of rotation in radians (pi/144 is default)
    graddist = zeros(2*pi/alpha,1);
    gradline = zeros(plc.vertex_distance+1,2*pi/alpha);
    gradlinerow = zeros(plc.vertex_distance+1,2*pi/alpha);
    gradlinecol = zeros(plc.vertex_distance+1,2*pi/alpha);
    parameterline = zeros(plc.vertex_distance+1,2*pi/alpha);
    i_coords = zeros(2*pi/alpha,4); % XC line vertex coordinates
    ix1 = zeros(2*pi/alpha,1);
    iy1 = zeros(2*pi/alpha,1);
    ix2 = zeros(2*pi/alpha,1);
    iy2 = zeros(2*pi/alpha,1);
    for i = 1:(2*pi/alpha)
        beta = i*alpha;
        plc.vxcenter = plc.max_gradient_col;
        plc.vycenter = plc.max_gradient_row;
        x1 = plc.vxcenter;
        x2 = plc.vxcenter;
        y1 = plc.vycenter + (plc.vertex_distance/2);
        y2 = plc.vycenter - (plc.vertex_distance/2);
        xcv = [x1,y1;...
              x2,y2];
        origin = [plc.vxcenter,plc.vycenter];
        ix1(i) = (((xcv(1,1)-origin(1,1))*cos((i-1)*beta)) - ...
            ((xcv(1,2)-origin(1,2))*sin((i-1)*beta)))+origin(1,1);
        iy1(i) = (((xcv(1,2)-origin(1,2))*cos((i-1)*beta)) + ...
            ((xcv(1,1)-origin(1,1))*sin((i-1)*beta)))+origin(1,2);
        ix2(i) = (((xcv(2,1)-origin(1,1))*cos((i-1)*beta)) - ...
            ((xcv(2,2)-origin(1,2))*sin((i-1)*beta)))+origin(1,1);
        iy2(i) = (((xcv(2,2)-origin(1,2))*cos((i-1)*beta)) + ...
            ((xcv(2,1)-origin(1,1))*sin((i-1)*beta)))+origin(1,2);
        i_coords = [ix1, iy1, ix2, iy2];

        graddist(i,1) = round(pdist([ix1(i),iy1(i);...
            ix2(i),iy2(i)],'euclidean'));

        [gradlinecol(:,i),gradlinerow(:,i),gradline(:,i)] = improfile(...
            plc.parameter_gradient,[ix1(i),ix2(i)],[iy1(i),iy2(i)],...
            graddist(i)+1,'bilinear'); % find the values along 
            % the line of cross section for the gradient
        [~,~,parameterline(:,i)] = improfile(...
            plc.gridded_data,[ix1(i),ix2(i)],[iy1(i),iy2(i)],...
            graddist(i)+1,'bilinear'); % find the values along 
            % the line of cross section
        gradline(:,i) = smooth(gradline(:,i),'moving');
    end
%%% END GRADIENT LINE ROTATION


%             gradlinesum(i) = peak2peak(plc.target_data(gradlinerow(:,i),gradlinecol(:,i)));
gradlinesum = peak2peak(gradline);
gradparamline = gradient(parameterline);
%             [~,columns] = find(gradlinesum == max(gradlinesum));
[~,columns] = find(gradparamline == max(abs(gradparamline)));
col = columns(1,1);
plc.gradient_cutline = gradline(:,col);
plc.gradline_col = gradlinecol(:,col);
plc.gradline_row = gradlinerow(:,col);

if isempty(plc.gradient_cutline) == 0
    gradient_vector = plc.parameter_gradient(:);
    gradient_vector = sort(gradient_vector(:),'descend','MissingPlacement','last');
    legitimate_value = false;
    k = 1;
    while legitimate_value == false
        sprintf('%s',k);
        maxgrad = gradient_vector(k,1); % identify the maximum gradient value
        [maxgradrow,maxgradcol] = find(plc.parameter_gradient == maxgrad);
            % find the location of the maxmimum gradient value
        plc.max_gradient_row = maxgradrow;
        plc.max_gradient_col = maxgradcol;
        plc.vxcenter = maxgradcol;
        plc.vycenter = maxgradrow;
        [gradline, gradlinecol, gradlinerow, graddist, ~] = plc.GradientLineRotation();
        [~,columns] = find(gradlinesum == max(gradlinesum));
        col = columns(1,1);
        plc.gradient_cutline = gradline(:,col);
        plc.gradline_col = gradlinecol(:,col);
        plc.gradline_row = gradlinerow(:,col);
        if isempty(plc.gradient_cutline) == 0
            if k > 30
                sprintf('%s','WARNING: ',num2str(k-1),...
                    ' cut lines had to be discarded for plcect ',...
                    plc.fname_components{plc.varmag_component},...
                    '. This could be due to a complicated ',...
                    'geometry, or a large cut line length.')
            end
            legitimate_value = true;
        else
            k = k + 1;                    
        end
    end
end
plc.grad_vertex_distance = graddist(col);
plc.vx1 = i_coords(col,1);
plc.vy1 = i_coords(col,2);
plc.vx2 = i_coords(col,3);
plc.vy2 = i_coords(col,4);

i = 1;