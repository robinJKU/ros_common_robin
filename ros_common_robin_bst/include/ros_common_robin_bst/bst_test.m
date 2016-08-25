%% bst test file
% Alexander Reiter, Institute of Robotics, JKU
% November 2015

clear all
close all
clc

%addpath('H:\DA\NURBS');

test_case = 5;

switch test_case
    case 1 % GM UE 2013S, Aufgabe 7: polynomiale Bezierkurve
        degree = 1;
        ctrl_pts = [0, 2, 3, 4, 6];
        n_ctrl_pts = size(ctrl_pts,2);
    case 2
        ctrl_pts = [0, 0.5, 5];
        degree = 2;
        knots = [zeros(1,degree+1), ones(1,degree+1)];
    case 3 % approximation
        degree = 7;
        t = 0:0.01:1;
        approx.par = t;
        approx.val = t.^2.*sqrt(1-t.^2);
        approx.der = zeros(size(approx.par));
        n_ctrl_pts = 50;
    case 4
        degree = 3;
        ctrl_pts = [0, 0.45, 1.4, 0.45, 0];
        n_ctrl_pts = size(ctrl_pts,2);
    case 5
        degree = 3;
        ctrl_pts = [1.2, 1.0, 0, -1, -1.2];
        n_ctrl_pts = size(ctrl_pts,2);
end
if exist('approx', 'var')
    spl = bst(degree, approx, n_ctrl_pts);
else
    if exist('knots', 'var')
        spl = bst(degree, ctrl_pts, knots);
    else
        spl = bst(degree, ctrl_pts);
    end
end
par = linspace(0,1,999);

figure()
n_plots = degree + 2;
for i = 1:n_plots
    subplot(n_plots,1,i)
    der = (i-1)*ones(size(par));
    plot(par, bst(spl, par, der));
    switch (i-1) 
        case 0
            title(sprintf('spline of degree %d', degree))
            ylabel('function');
        case 1
            ylabel('1-st derivative'); 
        case 2
            ylabel('2-nd derivative');
        case 3
            ylabel('3-rd derivative');
        otherwise
            ylabel(sprintf('%d-th derivative', (i-1))); 
    end
    grid on
end
xlabel('parameter')
% plot(linspace(min(t),max(t),CP), ctrl_pts(1,:), 'k+')
% max(abs(vals(1:end-1)-p'))
% subplot(3,1,2)
% der1 = diff(p)./diff(t(1:end-1));
% plot(t(1:end-2), der1)
% hold on
% grid on
% vals_der1 = bst(spline, t, ones(size(t)));
% plot(t, vals_der1, 'r');
% subplot(3,1,3)
% der2 = diff(der1)./diff(t(1:end-2));
% plot(t(1:end-3), der2)
% hold on
% grid on
% vals_der2 = bst(spline, t, 2*ones(size(t)));
% plot(t, vals_der2, 'r');

return
% 
%     case 2 % GM UE 2013S, Aufgabe 6
%         ctrl_pts = [-1, 0, 1;
%              1, 2, 1];
%         CP = size(ctrl_pts,2);
%         w = [ 1 1 1 ];
%         n = 2;
%         d = n;
%         m = size(ctrl_pts,2) + n + 1;
%         t = 0:0.01:1;
%         tk = [zeros(1,n+1), ones(1,n+1)];
%     case 3 % GM UE 2013S, Aufgabe 5
%         ctrl_pts = [-1, 0, 3;
%              3, 1, 2];
%         CP = size(ctrl_pts,2);
%         n = 2;
%         m = size(ctrl_pts,2) + n + 1;
%         w = ones(1,m-n-1);
%         d = n;
%         t = 0:0.01:1;
%         tk = [zeros(1,n+1), ones(1,n+1)];
%     case 4 % GM UE 2013S, Aufgabe 8a
%         ctrl_pts = [-1/70, 421/630, 3011/630, 489/70;
%              -19/70, 671/70, -1277/210, 51/70];
%         CP = size(ctrl_pts,2);
%         n = 2;
%         m = size(ctrl_pts,2) + n + 1;
%         w = ones(1,m-n-1);
%         d = n;
%         t = 0:0.01:1;
%         tk = [zeros(1,n+1), ones(1,n+1)];
%         %tk = linspace(0,1,2*n+2);
%     case 5 % GM UE 2013S, Aufgabe 8b
%         ctrl_pts = [0, -1/12, 67/18, 19/4, 7;
%              0, -47/12, 373/18, -185/12, 1];
%         n = 4;
%         CP = size(ctrl_pts,2);
%         m = size(ctrl_pts,2) + n + 1;
%         w = ones(1,m-n-1);
%         d = 4;
%         t = 0:0.01:1;
%         tk = [zeros(1,n+1), ones(1,n+1)];
%     case 6 % own test
%         n = 9;
%         ctrl_pts = linspace(0,1,n+1);
%         ctrl_pts = [ctrl_pts; 0, 0, (-1).^(1:n+1-4), 0, 0]; 
%         ctrl_pts(2,3) = 1000;
%         ctrl_pts(2,4) = -2000; 
%         ctrl_pts(2,5) = 1000;
%         CP = size(ctrl_pts,2);
%         m = size(ctrl_pts,2) + n + 1;
%         w = ones(1,m-n-1);
%         d = n;
%         t = 0:0.01:1;
%         tk = [zeros(1,n+1), ones(1,n+1)];
%     case 7
%         n = 6;
%         d = 6;
%         CP = 30;
%         ctrl_pts = [zeros(1,d), linspace(0,1,CP-2*d),ones(1,d)];
%         %c = [c;c];
%         m = CP+n+1;
%         w = ones(1,m-n-1);
%         tk = [zeros(1,d+1), 0.5+0.9*linspace(-.5,.5,m-2*d-2), ones(1,d+1)]-0.5;
%         t = (0:0.001:1)-0.5;
%     case 8
%         n = 2;
%         d = 2;
%         CP = 18;
%         ctrl_pts = [linspace(0,2*pi,CP); 0.5 + sin(linspace(0,2*pi,CP))];
%         ctrl_pts(2,10) = 4;
%         m = CP+n+1;
%         w = ones(1,m-n-1);
%         tk = [zeros(1,d+1), 0.5+0.8*linspace(-.5,.5,m-2*d-2), ones(1,d+1)];
%         t = 0:0.001:1;
%     case 9
%         n = 2;
%         d = 2;
%         ctrl_pts = [0, 1, 3, 4, 5;
%              0, 1, 2, 1, -1];
%         CP = size(ctrl_pts,2);
%         w = [1,1,1,1,1];
%         tk = [0,0,0,1/3,2/3,1,1,1];
%         m = length(tk);
%         t = 0:1e-3:1;
%     case 10
%         load('C:\Users\Alexander Reiter\Documents\MyJabberFiles\ak114689@jku.at\fCOP_Ident.mat');
%         n = SplineGrad_fCOP;
%         d = n;
%         ctrl_pts = ContPoints_fCOP;
%         CP = size(ctrl_pts,2);
%         w = ones(1,25);
%         tk = knots_fCOP;
%         m = length(tk);
%         t = min(tk):1e-4:max(tk);
% end
% tic
% %[p, N_all] = nurbs_new(n, d, w, c(1,:), tk, t );
% toc
% tic
% spline = bst(n, ctrl_pts(1,:), tk);
% toc
% %return
% figure;
% subplot(3,1,1)
% %plot(t(1:end), p, 'LineWidth', 2)
% hold on
% grid on
% tic
% vals = bst(spline, t);
toc
plot(t, vals, 'r');
plot(linspace(min(t),max(t),CP), ctrl_pts(1,:), 'k+')
max(abs(vals(1:end-1)-p'))
subplot(3,1,2)
der1 = diff(p)./diff(t(1:end-1));
plot(t(1:end-2), der1)
hold on
grid on
vals_der1 = bst(spline, t, ones(size(t)));
plot(t, vals_der1, 'r');
subplot(3,1,3)
der2 = diff(der1)./diff(t(1:end-2));
plot(t(1:end-3), der2)
hold on
grid on
vals_der2 = bst(spline, t, 2*ones(size(t)));
plot(t, vals_der2, 'r');


