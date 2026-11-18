function u = solve_robin_diffusion(un, face_diffusivity, source, dt, dy, transfer, ambient)
%SOLVE_ROBIN_DIFFUSION Conservative backward Euler with two Robin edges.
%   u_t = d/dy(D du/dy) + source; outward flux = transfer*(u-ambient).
%   un contains ALL physical nodes, including both gouge edges.
%   face_diffusivity is a scalar or ny-1 inter-node values [m^2/s].
%   source is a scalar or ny nodal values [field units/s].
%   transfer [m/s] and ambient are scalars or [south north] pairs.
%   transfer=0 means no flux; transfer=Inf imposes the ambient value.
%
% Half-cell boundary volumes (equivalently, eliminated centered ghost
% nodes) retain storage and sources at the edges. With finite transfer,
% trapezoidal storage changes exactly by sources minus outward fluxes.
% The returned vector has the same orientation as un.

validateattributes(un, {'numeric'}, {'real','finite','vector'});
validateattributes(dt, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(dy, {'numeric'}, {'real','finite','scalar','positive'});
ny = numel(un);
if ny < 3
    error('robin:GridSize', 'At least three physical nodes are required.');
end
face_diffusivity = expand_values(face_diffusivity, ny-1, 'face_diffusivity');
if any(face_diffusivity < 0)
    error('robin:Diffusivity', 'Face diffusivities must be nonnegative.');
end
source = expand_values(source, ny, 'source');
ambient = expand_values(ambient, 2, 'ambient');
validateattributes(transfer, {'numeric'}, {'real','vector','nonnegative','nonnan'});
if isscalar(transfer), transfer = repmat(transfer,2,1); end
if numel(transfer) ~= 2
    error('robin:TransferSize', 'transfer must be a scalar or [south north] pair.');
end
transfer = transfer(:);

r = dt .* face_diffusivity ./ dy^2;
lower = [0; -r];
upper = [-r; 0];
% Boundary control volumes have half the interior width.
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

% Thomas elimination. Positive storage and nonnegative diffusion/exchange
% produce a strictly diagonally dominant M-matrix; no pivot clipping.
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

function values = expand_values(values, n, name)
validateattributes(values, {'numeric'}, {'real','finite','vector'}, mfilename, name);
if isscalar(values), values = repmat(values,n,1); end
if numel(values) ~= n
    error('robin:InputSize', '%s must be scalar or contain %d values.', name, n);
end
values = values(:);
end
