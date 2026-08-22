clear; close all; clc;
load(fullfile('..', 'post', 'fault_slip_data.mat'));
load(fullfile('..', 'post', 'fault_slip_field.mat'));
addpath("/Users/yuhanwang/Project/Localiza/Manuscript/ListOfFigs/DrosteEffect-BrewerMap-3.2.5.0")
secs_per_year = 3600*24*365;
W_rsf     = pi.^2 .* a ./ f0.^2 .* rho_c ./ Lambda .* (alpha_hy + alpha_th) ./ history_V;
W_rsf_fit = sqrt(2/pi) .* history_V ./ (max(history_gammadot')');
Pa_to_Mpa = 1/1e6;
s_to_ms   = 1000;
alpha_hy = 6.71e-6; 

% =========================================================
% ABSOLUTE PANEL SIZE IN INCHES  (figure size is derived)
% =========================================================
panel_W_in = 4.20;         % data-area width  (inches)
panel_H_in = 1.40;         % data-area height (inches)  →  exact 3:1

left_in   = 0.75;          % y-label + tick labels
top_in    = 0.40;          % above first panel
gap_in    = 0.65;          % between panels (x-tick + x-label room)
bot_in    = 0.50;          % below last panel
cb_gap_in = 0.08;          % panel right edge → colorbar left edge
cb_W_in   = 0.15;          % colorbar bar width
right_in  = 0.40;          % after colorbar

% Derived figure size in inches
fig_W_in = left_in + panel_W_in + cb_gap_in + cb_W_in + right_in;
fig_H_in = top_in  + 3*panel_H_in + 2*gap_in + bot_in;

% Normalized positions (derived from inches above)
lft  = left_in    / fig_W_in;
pw   = panel_W_in / fig_W_in;
ph   = panel_H_in / fig_H_in;
cbg  = cb_gap_in  / fig_W_in;
cbw  = cb_W_in    / fig_W_in;
top  = top_in     / fig_H_in;
gap  = gap_in     / fig_H_in;

y1 = 1 - top - ph;
y2 = y1 - gap - ph;
y3 = y2 - gap - ph;

cb_x = lft + pw + cbg;
% =========================================================

Vc = pi^2 * a / (f0*(f0-2*a)) * rho_c / Lambda * (alpha_th + alpha_hy) / h;
gammadot0 = Vc/h;
wavelength = 2*pi*sqrt(a*rho_c / (f0*Lambda) * (alpha_th + alpha_hy) ./ (f0 + 2*a) ./ gammadot0);
wavenumber = 2*pi/wavelength;
omega = sqrt(abs((16 * pi^2 ./ wavelength.^4 * (alpha_th - alpha_hy).^2 - ...
    8 * pi^2 ./ wavelength.^2 .* (alpha_th + alpha_hy) .* f0 / a .* f0 .* Lambda / rho_c .* gammadot0 + ...
    (f0/(a)).^2 + (f0*Lambda / rho_c).^2 .* gammadot0.^2)))/2;
v_migration = omega/wavenumber;

index_local_delocal = (log10(history_V) > log10(Vc))';
d = diff([0 index_local_delocal 0]);
start_idx_all = find(d == 1);
end_idx_all   = find(d == -1);
index_local_delocal_x = index_local_delocal;
index_local_delocal_x(start_idx_all(1):end_idx_all(1)) = 0;
start_idx = start_idx_all(1:end);
end_idx   = end_idx_all(1:end);
idx_pairs = [start_idx(:) end_idx(:)];

for i_inx = 1:10
    index_start = idx_pairs(i_inx,1) - 50;
    index_end   = idx_pairs(i_inx,2) + 50;
    for i_inner_inx = index_start:index_end
        W_rsf_fit(i_inner_inx) = localization_thickness(yy, history_gammadot(i_inner_inx,:));
    end

    close all;
    fig = figure('Color','w', 'Renderer','painters');  % painters = true vector output
    set(fig, 'InvertHardcopy','off');

    % Set figure size in inches — this is what exportgraphics will use exactly
    set(fig, 'Units','inches', 'Position', [1, 1, fig_W_in, fig_H_in]);
    set(fig, 'PaperUnits','inches', ...
             'PaperSize',     [fig_W_in, fig_H_in], ...
             'PaperPosition', [0, 0, fig_W_in, fig_H_in]);

    % ---- Create axes ----
    ax1 = axes('Parent', fig);
    ax2 = axes('Parent', fig);
    ax3 = axes('Parent', fig);

    pos1 = [lft, y1, pw, ph];
    pos2 = [lft, y2, pw, ph];
    pos3 = [lft, y3, pw, ph];
    set(ax1, 'Position', pos1);
    set(ax2, 'Position', pos2);
    set(ax3, 'Position', pos3);

    steps   = 1:length(timev(index_start:index_end));
    %t_ms   = (timev(index_start:index_end) - timev(index_start)) * secs_per_year * s_to_ms;
    y_norm = yy / h;

    % ===================== (1) history_p =====================
    seg_p = history_p(index_start:index_end,:)' .* Pa_to_Mpa;
    imagesc(ax1, steps, y_norm, seg_p);
    set(ax1, 'YDir','normal', 'TickLabelInterpreter','latex', 'FontSize',15);
    xlabel(ax1, 'Steps', 'Interpreter','latex', 'FontSize',15);
    ylabel(ax1, '$y/h$',               'Interpreter','latex', 'FontSize',15);
    box(ax1,'on');
    pmax = max(seg_p(:));  pmin = min(seg_p(:));
    if pmax == 0, pmax = 1; end
    colormap(ax1, makeBuReMap(256, -20, pmax));
    clim(ax1, [-20, pmax]);
    cb1 = colorbar(ax1);
    set(ax1, 'Position', pos1);
    cb1.Units = 'normalized';
    cb1.Position = [cb_x, y1, cbw, ph];
    cb1.TickLabelInterpreter = 'latex';
    cb1.Title.String         = '$p \ [\mathrm{MPa}]$';
    cb1.Title.Interpreter    = 'latex';
    cb1.Title.FontSize       = 15;

    % ===================== (2) history_gammadot =====================
    imagesc(ax2, steps, y_norm, log10(history_gammadot(index_start:index_end,:)'));
    set(ax2, 'YDir','normal', 'TickLabelInterpreter','latex', 'FontSize',15);
    xlabel(ax2, 'Steps', 'Interpreter','latex', 'FontSize',15);
    ylabel(ax2, '$y/h$',               'Interpreter','latex', 'FontSize',15);
    box(ax2,'on');
    colormap(ax2, inferno(256));
    clim(ax2, [-4, 4]);
    cb2 = colorbar(ax2);
    set(ax2, 'Position', pos2);
    cb2.Units = 'normalized';
    cb2.Position = [cb_x, y2, cbw, ph];
    cb2.TickLabelInterpreter = 'latex';
    cb2.Title.String         = '$\log_{10}\dot{\gamma}$';
    cb2.Title.Interpreter    = 'latex';
    cb2.Title.FontSize       = 15;

    % ===================== (3) history_phi_p =====================
    imagesc(ax3, steps, y_norm, history_alphahy(index_start:index_end,:)'./alpha_hy);
    set(ax3, 'YDir','normal', 'TickLabelInterpreter','latex', 'FontSize',15);
    xlabel(ax3, 'Steps', 'Interpreter','latex', 'FontSize',15);
    ylabel(ax3, '$y/h$',               'Interpreter','latex', 'FontSize',15);
    box(ax3,'on');
    colormap(ax3, inferno(256));
    %clim(ax3, [0, 1]);
    cb3 = colorbar(ax3);
    set(ax3, 'Position', pos3);
    cb3.Units = 'normalized';
    cb3.Position = [cb_x, y3, cbw, ph];
    cb3.TickLabelInterpreter = 'latex';
    cb3.Title.String         = '$\alpha_{hy}/\alpha_{hy0}$';
    cb3.Title.Interpreter    = 'latex';
    cb3.Title.FontSize       = 15;

    % ===================== EXPORT =====================
    outpng = fullfile('..','post', strcat('wrsf_comparison', num2str(i_inx), '_clean.pdf'));
    exportgraphics(fig, outpng, 'BackgroundColor','white', ...
        'ContentType','vector');   % vector PDF — no Resolution arg, no rasterisation
end

% ---- Local functions ----------------------------------------
function [W_eff, xc, sigma_eff] = localization_thickness(x, gammadot)
    x = x(:); gammadot = gammadot(:);
    gammadot(gammadot < 0) = 0;
    dx = x(2) - x(1);
    V  = sum(gammadot) * dx;
    xc = sum(x .* gammadot) * dx / V;
    sigma_eff = sqrt(sum((x - xc).^2 .* gammadot) * dx / V);
    W_eff = 2 * sigma_eff;
end

function cmap = inferno(n)
    if nargin < 1, n = 256; end
    anchors = [
        0.001462 0.000466 0.013866; 0.046915 0.030324 0.150164
        0.142378 0.046242 0.308553; 0.258234 0.038571 0.406485
        0.416331 0.090834 0.432943; 0.578304 0.148039 0.404411
        0.735683 0.215906 0.329812; 0.865006 0.316822 0.226055
        0.954506 0.468744 0.099874; 0.987622 0.645320 0.039886
        0.964394 0.843848 0.273391; 0.988362 0.998364 0.644924];
    cmap = interp1(linspace(0,1,size(anchors,1)), anchors, linspace(0,1,n));
end

function BuReMap = makeBuReMap(n, cmin, cmax)
    if nargin < 1 || isempty(n),    n    = 1024; end
    if nargin < 2 || isempty(cmin), cmin = -20;  end
    if nargin < 3 || isempty(cmax), cmax = 100;  end
    if cmin >= cmax,  error('cmin must be < cmax.'); end
    if ~(cmin <= 0 && cmax >= 0), error('0 must lie between cmin and cmax.'); end
    nNeg = max(round(n * (-cmin)/(cmax-cmin)), 1);
    nPos = max(n - nNeg, 1);
    negColors = flipud(brewermap(nNeg,'Blues')); negColors(end,:) = [1 1 1];
    posColors = brewermap(nPos,'Reds');          posColors(1,:)   = [1 1 1];
    BuReMap = [negColors; posColors];
end