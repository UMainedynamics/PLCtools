classdef PLCobject
   properties
       microstructure
       parameter
       tempInterval
       stressStrainInterval
       kValue
       results % The target data for the parameter of interest is stored here
       results_gradient % The gradient of the target data for the parameter of interest is stored here
       datatable
       grains
       phases
       phasesnulled
       msResolution
       xgrid
       ygrid
       gridded_size
       total_size
       fname_components
       theta
       coords
       vx1
       vy1
       vx2
       vy2
       vxcenter
       vycenter
       xc_x
       xc_y
       cut_line_v1
       cut_line_v2
       cutlinelength
       cutline
       max_gradient
       max_gradient_row
       max_gradient_col
       vertex_distance
       grad_vertex_distance
       grad_x
       grad_y
       gradient_cutline
       gradline_row
       gradline_col
   end
   
   
   methods
       function obj = PLCobject(parameter, tempInterval, ...
               stressStrainInterval, varargin)
        % The PLCobject class is used to create objects which are used to
        % access and manipulate the data produced by the PLC toolbox. To
        % create a new instance of the PLCobject class, the following input
        % arguments are required:
        %
        % The parameter of interest ('Stress' 'Strain' 'Viscosity' or 
        % 'PowerDissipationDensity'), temperature interval, stress-strain 
        % interval, and (if using 'Stress' or 'Strain' as the parameter
        % of interest), the k-value.
        %
        % Optional parameters include 'resolution' , 'nullPhases' , and
        % 'boundaryFitness'. The 'resolution' parameter can be used to
        % change the resolution of the grid onto which the target data is
        % interpolated (default = 150; increase for higher resolution). The
        % 'nullPhases' parameter can be used to null the grains associated
        % with a particular phase number. This can be useful if you are
        % interested in calculating the results of only the matrix or only
        % grains of a particular phase. The 'boundaryFitness' parameter
        % refers to the aggressiveness or tightness of the fit of the
        % nulled zone when using the 'nullPhases' option. The default
        % boundaryFitness is 1, which is the most agressive / tightest fit,
        % but this can be reduced to a number closer to 0 if the boundary
        % resolution is too aggressive around complex geometries. For more
        % information, type in 'doc boundary' in the MATLAB Command Window
        % and look for the (s — Shrink factor) section.
        %
        % EXAMPLES (all examples use a temperature interval of 1 and a
        % stress-strain interval of 2):
        % 
        %   EXAMPLE 1: myObject = PLCobject('Viscosity',1,2);
        %   EXAMPLE 2: myObject = PLCobject('Stress',1,2,'kValue',10);
        %   EXAMPLE 3: myObject = PLCobject('Strain',1,2,'kValue',5,...
        %       'nullPhases',2,'boundaryFitness',0.85);
            expectedParameters = {'Stress','Strain','Viscosity', ...
                'PowerDissipationDensity'};
            defaultResolution = 150;
            defaultNullPhases = NaN;
            defaultBoundaryFitness = 1;
            p = inputParser;
                validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
                validKvalues = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 11);
                validFitness = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 1);
                addRequired(p,'parameter', ...
                    @(x) any(validatestring(x,expectedParameters)));
                addRequired(p,'tempInterval',validScalarPosNum);
                addRequired(p,'stressStrainInterval',validScalarPosNum)
                addParameter(p,'kValue', validKvalues);
                addParameter(p,'resolution',defaultResolution, ...
                    validScalarPosNum);
                addParameter(p,'nullPhases', defaultNullPhases, ...
                    validScalarPosNum);
                addParameter(p,'boundaryFitness', defaultBoundaryFitness, ...
                    validFitness);
            parse(p, parameter, tempInterval, ...
               stressStrainInterval, varargin{:});
            obj.parameter = p.Results.parameter;
            obj.tempInterval = p.Results.tempInterval;
            obj.stressStrainInterval = p.Results.stressStrainInterval;
            obj.kValue = p.Results.kValue;
            obj.msResolution = p.Results.resolution;
            obj.phasesnulled = p.Results.nullPhases;
            boundaryFitness = p.Results.boundaryFitness;
            
            % Load the data from the microstructure field
            obj = loadPLCdata(obj);
            
            % Retrieve grain and phase info from the microstructure field
            obj.grains = obj.microstructure.Grains;
            obj.phases = obj.microstructure.ElementPhases;
            obj.coords = obj.microstructure.MicroCoordinates / ...
                max(obj.microstructure.MicroCoordinates(:));
            data_map = [obj.coords, obj.results];
            obj.datatable = array2table(data_map, 'VariableNames', ...
                {'x','y',obj.parameter});
            F = scatteredInterpolant(obj.datatable.x,obj.datatable.y, ...
                obj.datatable.(obj.parameter));
            max_x = max(obj.datatable.x);
            max_y = max(obj.datatable.y);
            [obj.xgrid,obj.ygrid] = meshgrid(linspace(0,max_x,max_x*obj.msResolution), ...
                linspace(0,max_y,max_y*obj.msResolution));
            obj.results = F(obj.xgrid,obj.ygrid);
            obj.gridded_size = size(obj.results);
            obj.total_size = [max(obj.coords(:,1)), max(obj.coords(:,2))];

            if isnan(obj.phasesnulled) ~= 1
                [row_elem,~] = find(obj.phases == obj.phasesnulled);
                mst_mask = data_map(row_elem,:);
                mst_mask_table = array2table(mst_mask,...
                'VariableNames',{'x','y',obj.parameter});
                xv = mst_mask_table.x;
                yv = mst_mask_table.y;
                k = boundary(xv,yv,boundaryFitness);
                [in,on] = inpolygon(obj.xgrid,obj.ygrid,xv(k),yv(k));
                obj.results(in) = NaN;
                obj.results(on) = NaN;
            end
       end
       
       
       
       function obj = loadPLCdata(obj)
           % Load the microstructure data associated with the PLC results
           % in the present working directory.
           %
           % EXAMPLE: obj = loadPLCdata(obj);
            microstress = dir([pwd '/*_micro_stress_griddata.mat']);
            if isempty(microstress) == 1
                matfiles = dir([pwd '/*.mat']);
                matnames = {matfiles.name};
                namelength = cellfun('length',matnames);
                [~,col] = find(namelength == min(namelength));
                obj.microstructure = matnames{col};
            else
                run_name_delimiter = '_';
                obj.fname_components = strsplit(microstress.name, ...
                    run_name_delimiter);
                obj.microstructure = strcat(strjoin( ...
                    obj.fname_components(1:end-3),'_'),'.mat');
            end
            obj.microstructure = load(obj.microstructure);
            obj.microstructure = obj.microstructure.ms;
            obj.results = obj.microstructure.(strcat('Micro',obj.parameter));
            if strcmp(obj.parameter,'Vicsosity') || ...
                    strcmp(obj.parameter,'PowerDissipationDensity') == 1
                obj.results = obj.results{obj.tempInterval};
            else
                if strcmp(obj.parameter,'Stress') || ...
                        strcmp(obj.parameter,'Strain') == 1
                    obj.results = obj.results{obj.kValue,obj.tempInterval};
                else
                    error('k-value must be included if the parameter of interest is stress or strain')
                end
            end
            obj.results = obj.results(:,obj.stressStrainInterval);
       end
              function obj = CutLineGradient(obj)
           % Return the gradient cut line using manual vertices
            obj.gradient_cutline = smooth(improfile(obj.results_gradient, ...
                [obj.vx1,obj.vx2], [obj.vy1,obj.vy2], obj.vertex_distance+1, ...
                'bilinear'),'moving');  % find the values along the line
            obj.grad_vertex_distance = obj.vertex_distance;
       end

       
       
      function obj = SetCutLine(obj, cutline)
       % Sets the cut line vertices. Either provide a cut line length 
       % or provide vertices as an input argument as a 2x2 double 
       % in the format:s [x1,y1;x2,y2]
       %
       % EXAMPLE 1:
       %    cut_line_length = 0.2
       %    obj.SetCutLine(cut_line_length)
       %
       % EXAMPLE 2: 
       %    xcv = [0,0.5;...
       %          1,0.5]
       %    obj.SetCutLine(xcv)
       %
        if size(cutline) == 1
            obj = obj.CutLineMaxGradient(cutline);
        elseif size(cutline) == 2
            verts = cutline;
            obj.vx1 = obj.denormalize(verts(1,1));
            obj.vy1 = obj.denormalize(verts(1,2));
            obj.vx2 = obj.denormalize(verts(2,1));
            obj.vy2 = obj.denormalize(verts(2,2));
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
            xc_line = smooth(improfile(obj.results,...
                [obj.vx1,obj.vx2], [obj.vy1,obj.vy2], obj.vertex_distance+1,...
                'bilinear'),'moving');  % find the values along the line
            obj.cutline = xc_line;
       end
       
       
       
       function obj = CalculateGradient(obj)
       % Calculate the gradient of the parameter of interest as defined
       % by the obj.results data.
       %
       % EXAMPLE: obj(i) = obj(i).CalculateGradient
            [obj.grad_x,obj.grad_y] = gradient(obj.results,1);
            obj.results_gradient = sqrt((obj.grad_x).^2+(obj.grad_y).^2);       
            obj.max_gradient = max(obj.results_gradient(:));
       end
       
       
       function [gradline, gradlinerow, gradlinecol, graddist, parameterline, i_coords] = GradientLineRotation(obj)
       % Rotates the line of cross section about the centroid, which is
       % defined by the maximum gradient. Returns the rotated line, its
       % coordinate information, 
       %
       % EXAMPLE: [gradline, gradlinecol, gradlinerow, gradlinesum, graddist, i_coords] = obj.GradientLineRotation();
            alpha = pi/144; % angle of rotation in radians (pi/144 is default)
            graddist = zeros(2*pi/alpha,1);
            gradline = zeros(obj.cutlinelength+1,2*pi/alpha);
            gradlinerow = zeros(obj.cutlinelength+1,2*pi/alpha);
            gradlinecol = zeros(obj.cutlinelength+1,2*pi/alpha);
            parameterline = zeros(obj.cutlinelength+1,2*pi/alpha);
            i_coords = zeros(2*pi/alpha,4); % XC line vertex coordinates
            ix1 = zeros(2*pi/alpha,1);
            iy1 = zeros(2*pi/alpha,1);
            ix2 = zeros(2*pi/alpha,1);
            iy2 = zeros(2*pi/alpha,1);
                for i = 1:(2*pi/alpha)
                    beta = i*alpha;
                    obj.vxcenter = obj.max_gradient_col;
                    obj.vycenter = obj.max_gradient_row;
                    x1 = obj.vxcenter;
                    x2 = obj.vxcenter;
                    y1 = obj.vycenter + (obj.cutlinelength/2);
                    y2 = obj.vycenter - (obj.cutlinelength/2);
                    xcv = [x1,y1;...
                          x2,y2];
                    origin = [obj.vxcenter,obj.vycenter];
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
                        obj.results_gradient,[ix1(i),ix2(i)],[iy1(i),iy2(i)],...
                        graddist(i)+1,'bilinear'); % find the values along 
                        % the line of cross section for the gradient
                    [~,~,parameterline(:,i)] = improfile(...
                        obj.results,[ix1(i),ix2(i)],[iy1(i),iy2(i)],...
                        graddist(i)+1,'bilinear'); % find the values along 
                        % the line of cross section
                    gradline(:,i) = smooth(gradline(:,i),'moving');
                end
       end

           
       function obj = CutLineMaxGradient(obj, cutlinelength)
       % Solve for the cut line across the max gradient of the parameter of
       % interest.
            obj.cutlinelength = obj.denormalize(cutlinelength);
            [maxgradrow,maxgradcol] = find(obj.results_gradient == obj.max_gradient);
            obj.max_gradient_row = maxgradrow;
            obj.max_gradient_col = maxgradcol;
            obj.vxcenter = maxgradcol;
            obj.vycenter = maxgradrow;
            [gradline, gradlinecol, gradlinerow, graddist, parameterline, i_coords] = obj.GradientLineRotation();
            gradlinesum = peak2peak(gradline);
            gradparamline = gradient(parameterline);
            [~,columns] = find(gradparamline == max(abs(gradparamline)));
            col = columns(1,1);
            obj.gradient_cutline = gradline(:,col);
            obj.gradline_col = gradlinecol(:,col);
            obj.gradline_row = gradlinerow(:,col);
            if isempty(obj.gradient_cutline) == 0
                gradient_vector = obj.results_gradient(:);
                gradient_vector = sort(gradient_vector(:),'descend','MissingPlacement','last');
                legitimate_value = false;
                k = 1;
                while legitimate_value == false
                    sprintf('%s',k);
                    maxgrad = gradient_vector(k,1); % identify the maximum gradient value
                    [maxgradrow,maxgradcol] = find(obj.results_gradient == maxgrad);
                        % find the location of the maxmimum gradient value
                    obj.max_gradient_row = maxgradrow;
                    obj.max_gradient_col = maxgradcol;
                    obj.vxcenter = maxgradcol;
                    obj.vycenter = maxgradrow;
                    [gradline, gradlinecol, gradlinerow, graddist, ~] = obj.GradientLineRotation();
                    [~,columns] = find(gradlinesum == max(gradlinesum));
                    col = columns(1,1);
                    obj.gradient_cutline = gradline(:,col);
                    obj.gradline_col = gradlinecol(:,col);
                    obj.gradline_row = gradlinerow(:,col);
                    if isempty(obj.gradient_cutline) == 0
                        if k > 30
                            sprintf('%s','WARNING: ',num2str(k-1),...
                                ' cut lines had to be discarded for object ',...
                                obj.fname_components{obj.varmag_component},...
                                '. This could be due to a complicated ',...
                                'geometry, or a large cut line length.')
                        end
                        legitimate_value = true;
                    else
                        k = k + 1;                    
                    end
                end
            end
            obj.grad_vertex_distance = graddist(col);
            obj.vx1 = i_coords(col,1);
            obj.vy1 = i_coords(col,2);
            obj.vx2 = i_coords(col,3);
            obj.vy2 = i_coords(col,4);
       end
       
       
                   
       function max_less_min = MaxLessMin(obj, varargin)
       % Calculates the max-min of the parameter of interest across the line
       % of cross section. Unless a parameter of interest is provided as 
       % an input argument, the method defaults to using the parameter of
       % interest defined by the data provided in PLCvars.
       %
       % EXAMPLE: max_less_min(i,:) = obj.MaxLessMin(varmag_component);
            if size(varargin,2) == 1
                varmag_comp = varargin{1};
            else
                varmag_comp = obj.varmag_component;
            end
            max_less_min = zeros(1,2);
            max_param = max(obj.cutline);
            min_param = min(obj.cutline);
            param_x = char(obj.fname_components(varmag_comp));
            param_x = strrep(param_x, ',', '.');
            max_less_min(1,1) = str2double(param_x(end-2:end));
            max_less_min(1,2) = max_param - min_param;
       end
       
       
       function obj = RotateByTheta(obj, i, theta, xcv)
       % solve for new cut line vertices using the equations:
       %    x' = x*cos(theta) - y*sin(theta)
       %    y' = y*cos(theta) + x*sin(theta)
       % Expects xcv to be a 2x2 double with the initial x1,y1;x2,y2
       % EXAMPLE: xcv = [0,0.5;...
       %                1,0.5]
           if nargin == 4
            origin = [0.5,0.5];
            obj.theta = theta;
            obj.vx1 = obj.denormalize((((xcv(1,1)-origin(1,1))*cos((i-1)*theta)) - ...
                ((xcv(1,2)-origin(1,2))*sin((i-1)*theta)))+origin(1,1));
            obj.vy1 = obj.denormalize((((xcv(1,2)-origin(1,2))*cos((i-1)*theta)) + ...
                ((xcv(1,1)-origin(1,1))*sin((i-1)*theta)))+origin(1,2));
            obj.vx2 = obj.denormalize((((xcv(2,1)-origin(1,1))*cos((i-1)*theta)) - ...
                ((xcv(2,2)-origin(1,2))*sin((i-1)*theta)))+origin(1,1));
            obj.vy2 = obj.denormalize((((xcv(2,2)-origin(1,2))*cos((i-1)*theta)) + ...
                ((xcv(2,1)-origin(1,1))*sin((i-1)*theta)))+origin(1,2));
           else
               error('"i" , "theta" , and "xcv" must be included in the RotateByTheta method call')
           end
       end
       
       
       function var_norm = normalize(obj,var)
        % Convert native spatial units to 0-to-1 spatial units
        %
        % EXAMPLE:
        %     nx1 = obj.normalize(cut_line_length);
           var_norm = var / max(obj.gridded_size);
       end
       
       
       function var_norm = denormalize(obj,var)
        % Convert 0-to-1 spatial units to native spatial units
        %
        % EXAMPLE:
        %     dx1 = obj.denormalize(obj.vx1);
           var_norm = var * max(obj.gridded_size);
       end
      
      
       function plothandle = plot(obj, fighandle, plcData)
        % Creates a pseudocolor (checkerboard) plot of the results
        % associated with the object, such as obj.results or
        % obj.results_gradient. Requires a figure handle and the target
        % data as input parameters.
        %
        % EXAMPLE: 
        % myPlot = figure('units','normalized','Position',[0 0 1 1]);
        % obj.plot(myPlot,obj.results);
           plothandle = PLCplot(fighandle,obj,plcData);
       end
       
       
       function plothandle = subplot(obj, fighandle, m, n, p, plcData)
        % Creates a pseudocolor (checkerboard) subplot of the results
        % associated with the object, such as obj.results or
        % obj.results_gradient. Requires a figure handle, subplot
        % dimensions (see the MATLAB Help section on 'subplot' for more),
        % and the target data as input parameters.
        %
        % EXAMPLE: 
        % mySubplot = figure('units','normalized','Position',[0 0 1 1]);
        % parameterSubplot = obj.subplot(mySubplot,1,2,1,obj.results);
        % gradientSubplot = obj.subplot(mySubplot,1,2,2,obj.results_gradient);
           plothandle = PLCsubplot(fighandle,m,n,p,obj,plcData);
       end
   end
end
