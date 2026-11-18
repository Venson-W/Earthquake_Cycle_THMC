function stepper = make_p_T_stepper(fixed)
%MAKE_P_T_STEPPER Prepare the coupled diffusion step for repeated use.
%   STEPPER = MAKE_P_T_STEPPER(FIXED) validates the fixed parameters once.
%   FIXED fields: ny, alpha_th, rho_c, Lambda, dy, rect_m_frac, Pr, Er,
%   robin, temperature_ambient, pressure_ambient. Ambient values are scalar
%   or [south north] changes from the reference state.
%
%   [T,p] = STEPPER(Tn,pn,gamma,tauf,dt,M_alpha_hy,n_pl_term, ...
%                   use_plastic_term,decomposition)
%   uses the same arithmetic and source ordering as implicit_p_T_robin.
%   Dynamic arrays are checked for size, real finite values and admissible
%   diffusivity/mask values. Constant validation and transfer conversion
%   are kept outside the time-stepping loop.

validateattributes(fixed, {'struct'}, {'scalar'});
required = {'ny','alpha_th','rho_c','Lambda','dy','rect_m_frac','Pr','Er', ...
    'robin','temperature_ambient','pressure_ambient'};
for index = 1:numel(required)
    if ~isfield(fixed,required{index})
        error('robin:MissingSetting','Missing fixed.%s.',required{index});
    end
end
validateattributes(fixed.ny, {'numeric'}, {'real','finite','scalar','integer','>=',3});
validateattributes(fixed.alpha_th, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(fixed.rho_c, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(fixed.Lambda, {'numeric'}, {'real','finite','scalar'});
validateattributes(fixed.dy, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(fixed.rect_m_frac, {'numeric'}, {'real','finite','scalar'});
validateattributes(fixed.Pr, {'numeric'}, {'real','finite','scalar'});
validateattributes(fixed.Er, {'numeric'}, {'real','finite','scalar'});
[transfer_T,transfer_p] = robin_transfer_coefficients(fixed.robin,fixed.rho_c);
ambient_T = ambient_values(fixed.temperature_ambient,'temperature_ambient');
ambient_p = ambient_values(fixed.pressure_ambient,'pressure_ambient');
ny = fixed.ny;
alpha_th = fixed.alpha_th;
rho_c = fixed.rho_c;
Lambda = fixed.Lambda;
dy = fixed.dy;
rect_m_frac = fixed.rect_m_frac;
Pr = fixed.Pr;
Er = fixed.Er;
thermal_faces = repmat(alpha_th,ny-1,1);
stepper = @step;

    function [T,p] = step(Tn,pn,gamma,tauf,dt,M_alpha_hy,n_pl_term, ...
            use_plastic_term,decomposition)
        if ~isnumeric(Tn) || ~isnumeric(pn) || ~isvector(Tn) || ...
                numel(Tn) ~= ny || ~isequal(size(Tn),size(pn))
            error('robin:FieldShape','Tn and pn must have the same %d-node vector shape.',ny);
        end
        if ~isnumeric(tauf) || ~isscalar(tauf) || ~isnumeric(dt) || ~isscalar(dt)
            error('robin:DynamicScalar','tauf and dt must be numeric scalars.');
        end
        if ~isnumeric(M_alpha_hy) || ~isvector(M_alpha_hy)
            error('robin:FaceCount','M_alpha_hy must be a numeric vector.');
        end
        if numel(M_alpha_hy) == ny+1
            pressure_faces = M_alpha_hy(2:ny);
        elseif isscalar(M_alpha_hy)
            pressure_faces = repmat(M_alpha_hy,ny-1,1);
        elseif numel(M_alpha_hy) == ny-1
            pressure_faces = M_alpha_hy;
        else
            error('robin:FaceCount','M_alpha_hy must have 1, ny-1, or ny+1 values.');
        end
        gamma = nodes(gamma,'gamma');
        decomposition = nodes(decomposition,'decomposition');
        plastic = nodes(n_pl_term,'n_pl_term');
        plastic_switch = nodes(double(use_plastic_term),'use_plastic_term');
        dynamic_values = [Tn(:);pn(:);gamma;tauf;dt;M_alpha_hy(:); ...
            plastic;plastic_switch;decomposition];
        if ~isreal(dynamic_values) || any(~isfinite(dynamic_values))
            error('robin:NonfiniteState','All dynamic inputs must be real and finite.');
        end
        if dt <= 0
            error('robin:TimeStep','dt must be positive.');
        end
        if any(M_alpha_hy(:) < 0)
            error('robin:Diffusivity','Face diffusivities must be nonnegative.');
        end
        if any(plastic_switch ~= 0 & plastic_switch ~= 1)
            error('robin:PlasticSwitch','use_plastic_term must contain only zero or one.');
        end

        temperature_source = tauf .* gamma ./ rho_c - rect_m_frac .* Er .* decomposition;
        T = diffuse(Tn,thermal_faces,temperature_source,dt,transfer_T,ambient_T);
        pressure_source = Lambda .* (T(:)-Tn(:)) ./ dt + plastic_switch .* plastic ...
            + rect_m_frac .* Pr .* decomposition;
        p = diffuse(pn,pressure_faces(:),pressure_source,dt,transfer_p,ambient_p);
        if any(~isfinite(T(:))) || any(~isfinite(p(:)))
            error('robin:NonfiniteState','The coupled diffusion step produced a nonfinite field.');
        end
    end

    function values = nodes(values,name)
        if ~isnumeric(values) || ~isvector(values)
            error('robin:NodalSize','%s must be a numeric scalar or vector.',name);
        end
        if isscalar(values)
            values = repmat(values,ny,1);
        elseif numel(values) ~= ny
            error('robin:NodalSize','%s must be scalar or contain all %d nodal values.',name,ny);
        end
        values = values(:);
    end

    function u = diffuse(un,face_diffusivity,source,dt,transfer,ambient)
        % Keep expression ordering and Thomas operations identical to
        % solve_robin_diffusion; only invariant preparation is removed.
        r = dt .* face_diffusivity ./ dy^2;
        lower = [0; -r];
        upper = [-r; 0];
        upper(1) = 2 * upper(1);
        lower(end) = 2 * lower(end);
        diagonal = 1 - lower - upper;
        rhs = un(:) + dt .* source;
        for edge = 1:2
            if edge == 1, i = 1; else, i = ny; end
            if isinf(transfer(edge))
                lower(i) = 0;
                upper(i) = 0;
                diagonal(i) = 1;
                rhs(i) = ambient(edge);
            else
                exchange = 2 * dt * transfer(edge) / dy;
                diagonal(i) = diagonal(i) + exchange;
                rhs(i) = rhs(i) + exchange * ambient(edge);
            end
        end
        for i = 2:ny
            factor = lower(i) / diagonal(i-1);
            diagonal(i) = diagonal(i) - factor * upper(i-1);
            rhs(i) = rhs(i) - factor * rhs(i-1);
        end
        solution = zeros(ny,1);
        solution(end) = rhs(end) / diagonal(end);
        for i = ny-1:-1:1
            solution(i) = (rhs(i) - upper(i)*solution(i+1)) / diagonal(i);
        end
        u = reshape(solution,size(un));
    end
end

function ambient = ambient_values(ambient,name)
validateattributes(ambient, {'numeric'}, {'real','finite','vector'},mfilename,name);
if isscalar(ambient)
    ambient = repmat(ambient,2,1);
elseif numel(ambient) ~= 2
    error('robin:InputSize','%s must be scalar or a [south north] pair.',name);
end
ambient = ambient(:);
end
