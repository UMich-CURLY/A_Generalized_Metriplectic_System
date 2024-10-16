clear all;clc;close all;
%%
data = load("SO3_result4.mat");
plot_data = load("SO3_plot.mat");
% plot_data = plot_data.data;
%%
% subplot(2, 1, 1)
close all

fig_size = [16, 8];
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [5, 5, fig_size]);
subplot(1, 1, 1, 'Position', [0.1, 0.15, 0.85, 0.75])

hold on
box on
grid on
% plot(cost_sdp / 1e4 / 6000, '-o')
% plot(cost_lp / 1e4 / 6000, '-o')
plot(abs(data.cost_sdp(1:370))/ 1e4 ,  '-', "LineWidth", 2) %/ 1e4 / (length(data.data.y)) * 3)
plot(abs(data.cost_lp(1:370)) / 1e4 , '-.', "LineWidth", 2) % / 1e4 / (length(data.data.y)) * 3)
set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

xlabel("Iteration", "interpreter", "latex", "FontSize", 15)
ylabel("Loss", "interpreter", "latex", "FontSize", 15)
legend({"Metric Loss", "Entropy Loss"}, "interpreter", "latex")
xlim([0, 400])

set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("SO3_cvx", "-dpdf")
%%
% h = figure(2)
fig_size = [12, 12];
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [5, 5, fig_size]);
subplot(1, 1, 1, 'Position', [0.075, 0.075, 0.9, 0.9])
% subplot(1, 3, 1)
hold on
for k = 1:200 % length(plot_data.y)
    if false
        % plot(plot_data.y{k}(:, 1), plot_data.y{k}(:, 2), 'r-', 'LineWidth', 0.1)
    else
        plot3(plot_data.y{k}(1:5:end, 1), plot_data.y{k}(1:5:end, 2), plot_data.y{k}(1:5:end, 3), 'b-', 'LineWidth', 0.1)
    end
    plot3(plot_data.y{k}(1, 1), plot_data.y{k}(1, 2), plot_data.y{k}(1, 3), 'r.', "MarkerSize", 8)
end
% t = 0:0.01:pi*2;
% plot(cos(t), sin(2 * t), 'k-', "LineWidth", 2)
[X,Y,Z] = sphere(20);
surf(X, Y, Z, 'FaceAlpha', 0.1, "EdgeColor","none", "FaceColor", "b")
daspect([1, 1, 1])
% xlim([-2.5, 2.5])

xlabel("$p_x$", "Interpreter", "latex", "FontSize", 18)
ylabel("$p_y$", "Interpreter", "latex", "FontSize", 18)
zlabel("$p_z$", "Interpreter", "latex", "FontSize", 18)
view(35, 35)

box on
grid on


set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("SO3_traj", "-dpdf")
%%
% close all
figure(3)
set(gca, 'YScale', 'log')
% subplot(1, 3, 2)
hold on
for k = 1:20
    E = plot_data.S_fun(plot_data.y{k}')' + plot_data.H_fun(plot_data.y{k}')';
    % E(abs(E) < 1e-12) = 1e-12;
    % plot(plot_data.t{1}, -E)
    plot(plot_data.t{1}(2:end), E(2:end) - E(1:end-1) > 0)
end

box on
grid on

ylabel("$-E$", "Interpreter", "latex", "FontSize", 15)
xlabel("Time ($s$)", "Interpreter", "latex", "FontSize", 15)


% ylim([1e-6, 1e3])
% plot
% set(gca, 'XScale', 'log')