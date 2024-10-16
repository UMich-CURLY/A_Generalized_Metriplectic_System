clear all;clc;close all;
%%
data = load("linear_result4.mat");
plot_data = load("linear_plot.mat");
%%
% subplot(2, 1, 1)
close all

fig_size = [16, 8];
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [5, 5, fig_size]);
subplot(1, 1, 1, 'Position', [0.12, 0.15, 0.85, 0.75])

hold on
box on
grid on
% plot(cost_sdp / 1e4 / 6000, '-o')
% plot(cost_lp / 1e4 / 6000, '-o')
plot(abs(data.cost_sdp(1:694))/ 1e4 ,  '-', "LineWidth", 2) %/ 1e4 / (length(data.data.y)) * 3)
plot(abs(data.cost_lp(1:694))/ 1e4  , '-.', "LineWidth", 2) % / 1e4 / (length(data.data.y)) * 3)
set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

xlabel("Iteration", "interpreter", "latex", "FontSize", 15)
ylabel("Loss", "interpreter", "latex", "FontSize", 15)
legend({"Metric Loss", "Entropy Loss"}, "interpreter", "latex")
xlim([0, 710])

set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
% print("linear_cvx", "-dpdf")
%%
fig_size = [18, 12];
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [5, 5, fig_size]);
subplot(1, 1, 1, 'Position', [0.075, 0.12, 0.9, 0.87])
% subplot(1, 3, 1)
hold on
for k = 1:200 % length(plot_data.y)
    if false
        plot(plot_data.y{k}(:, 1), plot_data.y{k}(:, 2), 'r-', 'LineWidth', 0.1)
    else
        plot(plot_data.y{k}(1:5:end, 1), plot_data.y{k}(1:5:end, 2), 'b-', 'LineWidth', 0.8)
    end
    plot(plot_data.y{k}(1, 1), plot_data.y{k}(1, 2), 'r.', "MarkerSize", 8)
end
t = 0:0.01:pi*2;
plot(cos(t), sin(2 * t), 'k-', "LineWidth", 2)
daspect([1, 1, 1])
ylim([-2.5, 2.5] * 0.8)
xlim([-2.9, 2.9] * 1)

xlabel("$q$", "Interpreter", "latex", "FontSize", 18)
ylabel("$p$", "Interpreter", "latex", "FontSize", 18)

box on
grid on

set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("linear_traj", "-dpdf")
%%
% close all
figure(3)
set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% subplot(1, 3, 2)
hold on
for k = 1:20
    S = plot_data.S_fun(plot_data.y{k}')';
    H = plot_data.H_fun(plot_data.y{k}')';
    E = S + H;
    % E(abs(E) < 1e-12) = 1e-12;
    plot(plot_data.t{1}, -E)
    % plot(plot_data.t{1}, S)
    % plot(plot_data.t{1}, H)
end

box on
grid on

% ylim([1e-6, 1e3])
% plot
% set(gca, 'XScale', 'log')