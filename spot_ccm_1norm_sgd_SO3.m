clear;clc;close all;
%% get the conic constraints
n = 3;
kappa_metric = 2;
kappa_entropy = 6;

[sdp_prog_base, Gamma, vx, x, W_list] = getConicConstraints(n, kappa_metric);
[lp_prog_base , p, Jac, Jac_sym, xs]  = getPolyConstraints(n, kappa_entropy);
%% initialize
% x = msspoly('x', 2);
J = diag([0.1, 0.2, 0.15]);
H = 1 / 2 * x' * J^-1 * x;
S = - 1 / 2 * (x(1)^2 + x(2)^2 + x(3)^2 - 1)^2 - H;
% S = - 1 / 2 * (x(2)^2 + x(3)^2 + 4 * x(1)^4 - 4 * x(1)^2)^2 - H;
vx_ = monomials(x, 0:6);
[~, dm, M] = decomp([vx_; S]);
p_gt = full(M(end, :));

pnum = p_gt(:) + randn(size(p_gt(:))) * 1;
% pnum = randn(size(pnum));

data_raw = load("SO3_data.mat"); data_raw = data_raw.data;
% data = data_raw;
% data_raw.y = cell2mat(data_raw.y');
% data_raw.dy = cell2mat(data_raw.dy');

data = data_raw;

id = randperm(length(data_raw.y));
id = id(1:2000);

%%
proj_flag = false;
%%
cost_sdp = [];
cost_lp = [];

metric_log  = {};
entropy_log = {};
%%
for iter = 1:100000
    %%
    data = data_raw;
    
    id = randperm(length(data_raw.y));
    id = id(1:2000);

    data.y   = data_raw.y(id, :);
    data.dy  = data_raw.dy(id, :);
        
    dz   = data.dy;
    dHdx = data.y * J^-1;

    Omega_fun = skew_SO3(x);
    xop = data.y; 

    %% search the metric

    % min_K |z - K(., vv)|
    [sol_sdp, prog_sdp, sdp_vec] = searchMetric(pnum, dz, xop, dHdx, Gamma, Omega_fun, x, Jac, sdp_prog_base, proj_flag);
    % [sol_sdp, ~] = searchMetric(z, vv, Gamma, x, sdp_prog_base);

    %% search the entrSSSopy

    % min_p |z - Avv|
    Gamma_fun = sol_sdp.eval(Gamma);
    [sol_lp, prog_lp, lp_vec] = searchEntropy(p, dz, xop, dHdx, Gamma_fun, Omega_fun, x, Jac, lp_prog_base, proj_flag);
    % [sol_lp, ~] = searchEntropy(z, vv, p, Gamma_fun, x, Jac, lp_prog_base);
    %
    pnum = sol_lp.eval(p);
    %%
    cost_sdp(end+1) = double( sol_sdp.eval(sol_sdp.objective) );
    cost_lp(end+1) = double( sol_lp.eval(sol_lp.objective) );
    
    [~, ~, MM] = decomp([vx; Gamma_fun(:)]);

    metric_log{end+1} = MM;
    entropy_log{end+1} = double(pnum);
end
%%
close all

h = figure(1);

hold on
box on
grid on
cost_num = double(sol_lp.eval(lp_vec));

err_list = ["$\delta_x$", "$\delta_y$", "$\delta_z$"];

for k = 1:3
    subplot(3, 1, k)
    % plot(cost_num(k:3:end))
    plot(-sort(-abs(cost_num(k:3:end))))
    hold on
    plot(cost_num(k:3:end) * 0 + 1e-3, 'k--')
    if k == 1
        title("Element-wise regression loss")
    end

    if k == 3
        xlabel("Time(s)")
    end
    
    ylabel(err_list(k), "Interpreter", "latex", "FontSize", 14)
    set(gca, 'YScale', 'log')
    box on 
    grid on
end

set(h, 'Units','pixels');
% set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
% print("fig1", "-dpdf")

h = figure(2);

subplot(2, 1, 1)
hold on
box on
grid on
% plot(cost_sdp / 1e4 / 6000, '-o')
% plot(cost_lp / 1e4 / 6000, '-o')
plot(abs(cost_sdp / 1e4 / (length(data.y)) * 3), '-o')
plot(abs(cost_lp  / 1e4 / (length(data.y)) * 3), '-o')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlabel("Iteration")
ylabel("Loss")
legend({"LMI - metric", "LP - entropy"})

% subplot(2, 1, 2)
% hold on
% box on
% grid on
% cost_all = [cost_sdp,cost_lp];
% cost_all(1:2:end) = cost_sdp;
% cost_all(2:2:end) = cost_lp;
% plot((cost_all / 1e4 / 6000), '-o')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')




set(h, 'Units','pixels');
% set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
% print("fig2", "-dpdf")


%
% h =figure(3)
% Jw = data.w;% * data.J;
% plot3(Jw(:, 1), Jw(:, 2), Jw(:, 3), '.')
% box on
% grid on
% xlabel("$L_x$", "Interpreter","latex", "FontSize", 14)
% ylabel("$L_y$", "Interpreter","latex", "FontSize", 14)
% zlabel("$L_z$", "Interpreter","latex", "FontSize", 14)
% 
% set(h, 'Units','pixels');
% % set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
% % print("fig3", "-dpdf")