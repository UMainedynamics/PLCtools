classdef PLCobject
   properties
       parameter
       total_size
       fname_components
       varmag_component
       theta
       cutline
       vx1
       vy1
       vx2
       vy2
       vxcenter
       vycenter
       stress_file
       vertex_distance
       grad_vertex_distance
       parameter_index
       cut_line_v1
       cut_line_v2
       gradient
       gradient_sigma
       gradient_cutline
       matrix
       grains
   end
   methods
       function obj = PLCobject(PLCvars,varargin)
       % Expects PLCvars to be a cell containing the parameter of
       % interest, temperature interval, stress-strain interval, and
       % k-value. EXAMPLE: PLCvars = {'stress',1,1,10};
           obj.parameter = PLCvars{1};
           obj.stress_file = dir([pwd '/*' obj.parameter '_griddata.mat']);
           load(obj.stress_file.name);
           parameter_data = strcat('micro_',[obj.parameter]);
            if strcmp([obj.parameter],'visc') || strcmp([obj.parameter],'pdd') == 1
                obj.parameter_index = eval(strcat(parameter_data,'{',num2str(PLCvars{2}),...
                    ',',num2str(PLCvars{3}),'}'));
            else
                if size(PLCvars,2) > 3
                    obj.parameter_index = eval(strcat(parameter_data,'{',num2str(PLCvars{2}),...
                    ',',num2str(PLCvars{3}),',',num2str(PLCvars{4}),'}'));
                    obj.varmag_component = PLCvars{5};
                else
                    error('k-value must be included if "visc" or "pdd" is parameter of interest')
                end
            end
            obj.total_size = floor(size(obj.parameter_index));
            if size(varargin,2) >= 1
                run_name_delimiter = varargin{1};
            else
                run_name_delimiter = '_';
            end
            obj.fname_components = strsplit(obj.stress_file.name,run_name_delimiter);
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
            obj.vertex_distance = round(pdist([obj.vx1,obj.vy1;obj.vx2,obj.vy2],'euclidean'));           
            xc_line = smooth(improfile(obj.parameter_index,...
                [obj.vx1,obj.vx2], [obj.vy1,obj.vy2], obj.vertex_distance+1,...
                'bilinear'),'moving');  % find the values along the line
            obj.cutline = xc_line;
            obj.vertex_distance = obj.vertex_distance;
       end
       
       function max_less_min = MaxLessMin(obj, varargin)
       % calculate the max-min of the parameter of interest across the line
       % of cross section
            if size(varargin,2) == 1
                varmag_comp = varargin{1};
            else
                varmag_comp = obj.varmag_component;
            end
            max_less_min = zeros(1,2);
            max_param = max(obj.cutline);
            min_param = min(obj.cutline);
            param_x = char(obj.fname_components(varmag_comp));
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
       
       function obj = CalculateGradient(obj, kernel_value, varargin)
       % Calculate the gradient of the parameter of interest as defined
       % by the obj.parameter_index data. Expects a kernel value for the
       % calculation of the gradient as defined by Guanglei Xiong's
       % gaussgradient function.
       %
       % EXAMPLE 1: obj(i) = obj(i).CalculateGradient
       % EXAMPLE 2: obj(i) = obj(i).CalculateGradient([0.5,0.5])
       % EXAMPLE 3: obj(i) = obj(i).CalculateGradient([0.5,0.5],1e-02)
       % EXAMPLE 4: obj(i) = obj(i).CalculateGradient([0.5,0.5],1e-02, 1)

            if size(varargin,2) == 1 % If an argument is provided with cell
                % values corresponding to grains which should be masked
                % (presumably, all grains but matrix), subtract those
                % grains from obj.matrix
                msindex = varargin{1};
                msfname = strcat(strjoin(obj.fname_components(1:end-3),...
                    '_'),'.mat');
                msload = load(msfname);
                obj.grains = msload.ms.NumberGrains;
                GrainPolyLinesCell = msload.ms.GrainPolylines;
                GrainPolyLines = GrainPolyLinesCell{msindex{1,1},msindex{1,2}};
                xv = GrainPolyLines(:,1);
                yv = GrainPolyLines(:,2);
                k = boundary(xv,yv,1);
                [Y,X] = find(obj.parameter_index);
                [in,on] = inpolygon(X,Y,xv(k),yv(k));
                obj.matrix = obj.parameter_index;
                obj.matrix(in) = NaN;
                obj.matrix(on) = NaN;
                obj.matrix(obj.matrix==Inf) = 0;
                obj.matrix(obj.matrix==-Inf) = 0;
                
%                 threshold = varargin{1}/100;
%                 X = obj.parameter_index/max(obj.parameter_index(:));
%                 X(isnan(X)) = 0;
%                 X(X==Inf) = 0;
%                 X(X==-Inf) = 0;
%                 W = min(X(X>0));
%                 X((X-W)>threshold) = 0;
%                 maskedImage = X;
%                 
%                 % Determine which values fall outside of the grain boundary
%                 obj.matrix = obj.parameter_index;
%                 obj.matrix(maskedImage>0) = NaN;
%                 if sum(~isnan(obj.matrix)) < 1
%                     error('The current threshold value resulted in a null grain matrix')
%                 end
            end

            % Find the gradient of the matrix
            [gx,gy] = gaussgradient(obj.matrix,kernel_value);
            obj.gradient = abs(gx+gy);
       end
       
       function [gradline, gradlinesum, graddist, i_coords] = GradientLine(obj, cutlinelength)
          if nargin == 2
            cutlinelength = cutlinelength * (obj.total_size(1,2)-1);
            alpha = pi/144; % angle of rotation in radians
            graddist = zeros(2*pi/alpha,1);
            gradline = zeros(cutlinelength+1,2*pi/alpha);
            gradlinesum = zeros(1,2*pi/alpha);
            i_coords = zeros(2*pi/alpha,4);
            ix1 = zeros(2*pi/alpha,1);
            iy1 = zeros(2*pi/alpha,1);
            ix2 = zeros(2*pi/alpha,1);
            iy2 = zeros(2*pi/alpha,1);
                for i = 1:(2*pi/alpha)
                    beta = i*alpha;
                    x1 = (obj.vxcenter);
                    x2 = (obj.vxcenter);
                    y1 = (obj.vycenter + (cutlinelength/2));
                    y2 = (obj.vycenter - (cutlinelength/2));
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
                    graddist(i,1) = round(pdist([ix1(i),iy1(i);...
                        ix2(i),iy2(i)],'euclidean'));
                    gradline(:,i) = smooth(improfile(obj.gradient,...
                        [ix1(i),ix2(i)],[iy1(i),iy2(i)],graddist(i)+1,...
                        'bilinear'),'moving'); % find the values along the
                        % line of cross section
                    i_coords = [ix1, iy1, ix2, iy2];
                    if any(isnan(gradline(:,i))) == 1
                        gradline(:,i) = NaN;
                    end
                    gradlinesum(i) = sum(gradline(:,i));
                end
           else
               error('"linelength" must be included in the CutLineMaxGradient method call')
          end
       end
           
       function obj = CutLineMaxGradient(obj, cutlinelength)
       % Solve for the cut line across the max gradient of the parameter of
       % interest.
                     
           if nargin == 2
            param_vec = obj.gradient(:);
            param_vec = sort(param_vec(:),'descend','MissingPlacement','last');
            legitimate_value = false;
            k = 1;
            while legitimate_value == false
                sprintf('%s',k);
                maxgrad = param_vec(k,1);
                [maxgradrow,maxgradcol] = find(obj.gradient == maxgrad);
                obj.vxcenter = maxgradcol;
                obj.vycenter = maxgradrow;
                [gradline, gradlinesum, graddist, i_coords] = obj.GradientLine(cutlinelength);
                [~,col] = find(gradlinesum == max(gradlinesum));
                obj.gradient_cutline = gradline(:,col);
                if isempty(obj.gradient_cutline) == 0
                    if k > 10
                        sprintf('%s','WARNING: ',num2str(k-1),...
                            ' cut lines had to be discarded for object ',...
                            obj.fname_components{obj.varmag_component},...
                            '. Consider using a larger threshold value.')
                    end
                    legitimate_value = true;
                else
                    k = k + 1;                    
                end                
            end
            obj.grad_vertex_distance = graddist(col);
            obj.vx1 = i_coords(col,1);
            obj.vy1 = i_coords(col,2);
            obj.vx2 = i_coords(col,3);
            obj.vy2 = i_coords(col,4);
           else
               error('"cutlinelength" must be included in the CutLineMaxGradient method call')
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
   end
end
