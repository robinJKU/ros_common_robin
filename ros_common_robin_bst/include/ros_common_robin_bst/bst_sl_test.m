clear all
close all
clc

%% to be changed by user
% list of derivatives to be computed
max_der = 3;
SHOW_PLOTS = 1;
SHOW_ERROR = 1;

% init spline
degree = 6;
N_ctrl_pts = 30;
ctrl_pts = [zeros(1,degree), linspace(0,1,N_ctrl_pts-2*degree),ones(1,degree)];
t_vec = 0:0.001:1;
spline = bst(degree, ctrl_pts);

% Simulink
T = max(t_vec);
Ts = min(diff(t_vec));
der_mask = 0:max_der;
return
%%
sim('bst_sl_model')

if SHOW_PLOTS
    figure;
    for i = der_mask
        subplot(max_der+1,1,i+1)
        vals = bst(spline, t_vec, i*ones(1,length(t_vec)));
        plot(t_vec, simout_vals(:,i+1));
        ylabel(sprintf('Spline derivative %d', i));
        grid on
        hold all
        plot(t_vec, vals);
    end
    xlabel('t')
end

if SHOW_ERROR
    figure;
    for i = der_mask
        subplot(max_der+1,1,i+1)
        vals = bst(spline, t_vec, i*ones(1,length(t_vec)));
        plot(t_vec, vals - simout_vals(:,i+1));
        ylabel(sprintf('Error derivative %d', i, 'r'));
    end
    xlabel('t')
end

