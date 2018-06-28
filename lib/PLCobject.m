classdef PLCobject
   properties
       parameter
       tempInterval
       stressStrainInterval
       kValue
       microstructure
       total_size
       real_size
       varmag_component
       parameter_file
       theta
       target_data
       cut_line_v1
       cut_line_v2
       cutlinelength
       cutline
       vx1
       vy1
       vx2
       vy2
       vxcenter
       vycenter
       vertex_distance
       grad_vertex_distance
       parameter_gradient
       grad_x
       grad_y
       max_gradient
       max_gradient_row
       max_gradient_col
       gradient_sigma
       gradient_cutline
       gradline_row
       gradline_col
       gradline_coords_round
       matrix
       grains
       graincount
       phases
       phasescount
       phasesnulled
       null_coords
       nullspace
   end
   
   
   methods
       function obj = PLCobject(parameter, tempInterval, ...
               stressStrainInterval, varargin)
        % Requires as input arguments the parameter of interest 
        % ('stress' 'strain' 'visc' or 'pdd'), temperature interval, 
        % stress-strain interval, and k-value.
        %
        % EXAMPLE : loadPLCdata('parameter','stress', ...
        % 'tempInterval',1,'stressStrainInterval',1, ...
        % 'kValue',10);
            expectedParameters = {'stress','strain','visc','pdd'};
            
            p = inputParser;
                validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
                validKvalues = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 11);
%                 p.Results.FunctionName = 'loadPLCdata';
                addRequired(p,'parameter', ...
                    @(x) any(validatestring(x,expectedParameters)));
                addRequired(p,'tempInterval',validScalarPosNum);
                addRequired(p,'stressStrainInterval',validScalarPosNum)
                addParameter(p,'kValue', validKvalues);
            parse(p, parameter, tempInterval, ...
               stressStrainInterval, varargin{:});
            obj.parameter = p.Results.parameter;
            obj.tempInterval = p.Results.tempInterval;
            obj.stressStrainInterval = p.Results.stressStrainInterval;
            obj.kValue = p.Results.stressStrainInterval;
            % Load the data from the microstructure field
            loadPLCdata(obj);
            
            % Retrieve grain and phase info from the microstructure field
            obj.grains = obj.microstructure.Grains;
            obj.phases = obj.microstructure.GrainPhases;
            obj.graincount = size((obj.grains),2);
            % count the number of individual grains
            obj.phasescount = size(unique(obj.phases),1);
            % count distinct phases are in the list of grains
            obj.target_data = flipud(obj.target_data);
            obj.total_size = floor(size(obj.target_data));
       end
       
       
       function obj = loadPLCdata(obj)
            microstress = dir([pwd '/*_micro_stress_griddata.mat']);
            run_name_delimiter = '_';
            fname_components = strsplit(microstress.name, ...
                run_name_delimiter);
            obj.microstructure = strcat(strjoin( ...
                fname_components(1:end-3),'_'),'.mat');
            obj.microstructure = load(obj.microstructure);
            obj.microstructure = obj.microstructure.ms;

            parameter_data = strcat('micro_',obj.parameter);  
            if strcmp(obj.parameter,'visc') || strcmp(obj.parameter,'pdd') == 1
                obj.target_data = eval(strcat(parameter_data,'{', ...
                    num2str(obj.tempInterval),',', ...
                    num2str(obj.stressStrainInterval),'}'));
            else
                if strcmp(obj.parameter,'stress') || strcmp(obj.parameter,'strain') == 1
                    obj.target_data = eval(strcat(parameter_data,'{',...
                        num2str(obj.tempInterval),',', ...
                        num2str(obj.stressStrainInterval),',', ...
                        num2str(obj.kValue),'}'));
                else
                    error('k-value must be included if the parameter of interest is stress or strain')
                end
            end
        end
       
       
       function obj = CalculateGradient(obj, cutlinelength, varargin)
       % Calculate the gradient of the parameter of interest as defined
       % by the obj.target_data data. Expects a kernel value for the
       % calculation of the gradient as defined by Guanglei Xiong's
       % gaussgradient function. Adding another argument allows masking of
       % grains by phase number.
       %
       % EXAMPLE 1: obj(i) = obj(i).CalculateGradient
       % EXAMPLE 2: obj(i) = obj(i).CalculateGradient(1)
       % ---> this example explicitly defines the kernel value of 1
       % EXAMPLE 3: obj(i) = obj(i).CalculateGradient(1,2)
       % ---> this example masks the grains with a phase number of 2
            obj.cutlinelength = cutlinelength * (obj.total_size(1,2)-1);
            [obj.grad_x,obj.grad_y] = gradient(obj.target_data,1);
            obj.parameter_gradient = sqrt((obj.grad_x).^2+(obj.grad_y).^2);
            
            if size(varargin,2) == 1
                obj.phasesnulled = varargin{1};
                obj = obj.MaskPhases;
            end
            
            obj.max_gradient = max(obj.parameter_gradient(:));
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
                        obj.parameter_gradient,[ix1(i),ix2(i)],[iy1(i),iy2(i)],...
                        graddist(i)+1,'bilinear'); % find the values along 
                        % the line of cross section for the gradient
                    [~,~,parameterline(:,i)] = improfile(...
                        obj.target_data,[ix1(i),ix2(i)],[iy1(i),iy2(i)],...
                        graddist(i)+1,'bilinear'); % find the values along 
                        % the line of cross section
                    gradline(:,i) = smooth(gradline(:,i),'moving');
%                     gradlinesum(i) = peak2peak(obj.target_data(gradlinerow(:,i),gradlinecol(:,i)));
%                     gradlinesum(i) = peak2peak(gradline(:,i));
                end
       end

           
       function obj = CutLineMaxGradient(obj)
       % Solve for the cut line across the max gradient of the parameter of
       % interest.      
            [maxgradrow,maxgradcol] = find(obj.parameter_gradient == obj.max_gradient);
            obj.max_gradient_row = maxgradrow;
            obj.max_gradient_col = maxgradcol;
            obj.vxcenter = maxgradcol;
            obj.vycenter = maxgradrow;
            [gradline, gradlinecol, gradlinerow, graddist, parameterline, i_coords] = obj.GradientLineRotation();
%             gradlinesum(i) = peak2peak(obj.target_data(gradlinerow(:,i),gradlinecol(:,i)));
            gradlinesum = peak2peak(gradline);
            gradparamline = gradient(parameterline);
%             [~,columns] = find(gradlinesum == max(gradlinesum));
            [~,columns] = find(gradparamline == max(abs(gradparamline)));
            col = columns(1,1);
            obj.gradient_cutline = gradline(:,col);
            obj.gradline_col = gradlinecol(:,col);
            obj.gradline_row = gradlinerow(:,col);
                
            if isempty(obj.gradient_cutline) == 0
                gradient_vector = obj.parameter_gradient(:);
                gradient_vector = sort(gradient_vector(:),'descend','MissingPlacement','last');
                legitimate_value = false;
                k = 1;
                while legitimate_value == false
                    sprintf('%s',k);
                    maxgrad = gradient_vector(k,1); % identify the maximum gradient value
                    [maxgradrow,maxgradcol] = find(obj.parameter_gradient == maxgrad);
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
       
       
       function obj = MaskPhases(obj)
       % EXAMPLE: obj(i) = obj(i).MaskPhases(2)
       % ---> this example masks the grains with a phase number of 2
       
       % If an argument is provided with cell
       % values corresponding to phases, subtract those phases
       % from the area in which the gradient is calculated
            coords = ms.MicroCoordinates/max(ms.MicroCoordinates(:));
            mst =[coords, ms.MicroStress{10, 1}(:,end)];
            mst_table = array2table(mst, ...
                'VariableNames',{'x','y','stress'});
            [grainRow,~] = find(obj.phases==obj.phasesnulled);
            for i = grainRow'
                target_grain = [obj.grains{1,i}];
                xv = target_grain(:,1);
                yv = target_grain(:,2);
                k = boundary(xv,yv,0.95);
                [Y,X] = find(obj.target_data);
                [in,on] = inpolygon(X,Y,xv(k),yv(k));
                %nullpolygon = [xv(k),yv(k)];
                matrixspace = obj.target_data;
                matrixspace(in) = NaN;
                matrixspace(on) = NaN;
                matrixspace = flipud(matrixspace);
                obj.matrix = obj.target_data;
                obj.matrix(isnan(matrixspace)) = NaN;
            end
            obj.nullspace = isnan(obj.matrix)==1;
            [rownull,colnull] = find(obj.nullspace==1);
            obj.null_coords = [rownull,colnull];

            obj.matrix(obj.matrix==Inf) = 0;
            obj.matrix(obj.matrix==-Inf) = 0;
            obj.parameter_gradient(isnan(obj.matrix)==1) = NaN;
            gradline_row_round = round(obj.gradline_row);
            gradline_col_round = round(obj.gradline_col);
            obj.gradline_coords_round = [gradline_row_round,gradline_col_round];
       end

       
      function obj = SetCutLine(obj, varargin)
       % Either use the object's vertices as defined in another method
       % (such as CutLineMaxGradient) or provide vertices as an input
       % argument as a 2x2 double with the initial x1,y1;x2,y2
       %
       % EXAMPLE 1:
       %    obj.SetCutLine
       %
       % EXAMPLE 2: 
       %    xcv = [0,0.5;...
       %          1,0.5]
       %    obj.SetCutLine(xcv)
       %
           switch nargin           
               case 2
                verts = varargin{1};
                obj.vx1 = obj.Denormalize(verts(1,1),'x');
                obj.vy1 = obj.Denormalize(verts(1,2),'y');
                obj.vx2 = obj.Denormalize(verts(2,1),'x');
                obj.vy1 = obj.Denormalize(verts(2,2),'y');                  
           end
            obj.cut_line_v1 = [obj.vx1,obj.vy1]; % cartesian coordinates for
                % the 1st vertex of the cut line / line of cross section
            obj.cut_line_v2 = [obj.vx2,obj.vy2]; % cartesian coordinates for
                % the 2nd vertex of the cut line / line of cross section
            obj.vertex_distance = round(pdist([obj.vx1,obj.vy1;...
                obj.vx2,obj.vy2],'euclidean'));           
            xc_line = smooth(improfile(obj.target_data,...
                [obj.vx1,obj.vx2], [obj.vy1,obj.vy2], obj.vertex_distance+1,...
                'bilinear'),'moving');  % find the values along the line
            obj.cutline = xc_line;
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
            obj.vx1 = obj.Denormalize((((xcv(1,1)-origin(1,1))*cos((i-1)*theta)) - ...
                ((xcv(1,2)-origin(1,2))*sin((i-1)*theta)))+origin(1,1),'x');
            obj.vy1 = obj.Denormalize((((xcv(1,2)-origin(1,2))*cos((i-1)*theta)) + ...
                ((xcv(1,1)-origin(1,1))*sin((i-1)*theta)))+origin(1,2),'y');
            obj.vx2 = obj.Denormalize((((xcv(2,1)-origin(1,1))*cos((i-1)*theta)) - ...
                ((xcv(2,2)-origin(1,2))*sin((i-1)*theta)))+origin(1,1),'x');
            obj.vy2 = obj.Denormalize((((xcv(2,2)-origin(1,2))*cos((i-1)*theta)) + ...
                ((xcv(2,1)-origin(1,1))*sin((i-1)*theta)))+origin(1,2),'y');
           else
               error('"i" , "theta" , and "xcv" must be included in the RotateByTheta method call')
           end
       end
       
       
       function var_norm = Normalize(obj,var,dimension)
       % Convert native spatial units to 0-to-1 spatial units
       %
       % EXAMPLE:
       %     nx1 = obj.Normalize(obj.vx1,'x');
       %     ny1 = obj.Normalize(obj.vy1,'y');
       %     nx2 = obj.Normalize(obj.vx2,'x');
       %     ny2 = obj.Normalize(obj.vy2,'y');

           if strcmp(dimension,'X') == 1 || strcmp(dimension,'x') == 1
               var_norm = var / (obj.total_size(1,1)-1);
           elseif strcmp(dimension,'Y') == 1 || strcmp(dimension,'y') == 1
               var_norm = var / (obj.total_size(1,2)-1);
           else
               error('the dimension must be provided as X, x, Y, or y')
           end
       end
       
       
      function var_norm = Denormalize(obj,var,dimension)
      % Convert 0-to-1 spatial units to native spatial units
      %
      % EXAMPLE:
      %     dx1 = obj.Denormalize(obj.vx1,'x');
      %     dy1 = obj.Denormalize(obj.vy1,'y');
      %     dx2 = obj.Denormalize(obj.vx2,'x');
      %     dy2 = obj.Denormalize(obj.vy2,'y');

           if strcmp(dimension,'X') == 1 || strcmp(dimension,'x') == 1
               var_norm = var * (obj.total_size(1,1)-1);
           elseif strcmp(dimension,'Y') == 1 || strcmp(dimension,'y') == 1
               var_norm = var * (obj.total_size(1,2)-1);
           else
               error('the dimension must be provided as X, x, Y, or y')
           end
      end
      
      
      function obj = MaskByThreshold(threshold_value)
        threshold = threshold_value/100;
        X = obj.target_data/max(obj.target_data(:));
        X(isnan(X)) = 0;
        X(X==Inf) = 0;
        X(X==-Inf) = 0;
        W = min(X(X>0));
        X((X-W)>threshold) = 0;
        maskedImage = X;
        % Determine which values fall outside of the grain boundary
        obj.matrix = obj.target_data;
        obj.matrix(maskedImage>0) = NaN;
        if sum(~isnan(obj.matrix)) < 1
            error('The current threshold value resulted in a null grain matrix')
        end
      end
      
      
%        function obj = GradientLineTranslation(obj)
%             alpha = pi/36; % angle of rotation in radians (pi/144 is default)
%             graddist = zeros(2*pi/alpha,1);
%             gradline = zeros(obj.cutlinelength+1,2*pi/alpha);
%             gradlinesum = zeros(1,2*pi/alpha);
%             i_coords = zeros(2*pi/alpha,4);
%             ix1 = zeros(2*pi/alpha,1);
%             iy1 = zeros(2*pi/alpha,1);
%             ix2 = zeros(2*pi/alpha,1);
%             iy2 = zeros(2*pi/alpha,1);
%                 for i = 1:(2*pi/alpha)
%                     beta = i*alpha;
%                     x1 = obj.vxcenter;
%                     x2 = obj.vxcenter;
%                     y1 = obj.vycenter + (obj.cutlinelength/2);
%                     y2 = obj.vycenter - (obj.cutlinelength/2);
%                     xcv = [x1,y1;...
%                           x2,y2];
%                     origin = [obj.vxcenter,obj.vycenter];
%                     ix1(i) = (((xcv(1,1)-origin(1,1))*cos((i-1)*beta)) - ...
%                         ((xcv(1,2)-origin(1,2))*sin((i-1)*beta)))+origin(1,1);
%                     iy1(i) = (((xcv(1,2)-origin(1,2))*cos((i-1)*beta)) + ...
%                         ((xcv(1,1)-origin(1,1))*sin((i-1)*beta)))+origin(1,2);
%                     ix2(i) = (((xcv(2,1)-origin(1,1))*cos((i-1)*beta)) - ...
%                         ((xcv(2,2)-origin(1,2))*sin((i-1)*beta)))+origin(1,1);
%                     iy2(i) = (((xcv(2,2)-origin(1,2))*cos((i-1)*beta)) + ...
%                         ((xcv(2,1)-origin(1,1))*sin((i-1)*beta)))+origin(1,2);
%                     i_coords = [ix1, iy1, ix2, iy2];
%                     
%                     graddist(i,1) = round(pdist([ix1(i),iy1(i);...
%                         ix2(i),iy2(i)],'euclidean'));
%                     
%                     [gradlinecol,gradlinerow,gradline(:,i)] = improfile(...
%                         obj.parameter_gradient,[ix1(i),ix2(i)],[iy1(i),iy2(i)],...
%                         graddist(i)+1,'bilinear'); % find the values along 
%                         % the line of cross section
%                     gradline(:,i) = smooth(gradline(:,i),'moving');
%                     gradlinesum(i) = sum(gradline(:,i));
%                 end
%        end

      
   end
end
