clearvars;
close all;
clc;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Arial Cyr');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Arial Cyr');
LW = 'LineWidth';
lw = 2;

% u_tt = u_xx + f;

%% Space initialization
Lx = 1;
nx = 300;
x = linspace(0, Lx, nx);
dx = x(2) - x(1);

%% Time initialization
T = 1;
c = 1;
CFL = 1;
dt = CFL * dx / c;
t_arr = 0:dt:T;
nt = length(t_arr);

%% Solution initialization
u_num_sol = zeros(nt, nx);
u_analyt_sol = zeros(nt, nx);
u_n = zeros(nx, 1);
u_nm1 = u_n; 
u_np1 = u_n;

%% Initial conditions
u_n = 0.5*sin(2*pi * x);
u_np1(:) = u_n(:);

%% Solution
t = 0;
figure('Color','White');
info_dat = zeros(nx * nt, 3);
info_dat_analyt = zeros(nx * nt, 3);
ch = 0;
while t < T
    % Reflecting
    u_n([1 end]) = 0;
    
    % Absorbing
%     u_np1(1) = u_n(2)  +((CFL - 1) / (CFL + 1)) * (u_np1(2) - u_n(1));
%     u_np1(end) = u_n(end - 1) + ((CFL - 1) / (CFL + 1)) * (u_np1(end - 1) - u_n(end));
        
    analyt_solve = 0.5*sin(2*pi*x) * cos(2*pi*t);
    
    info_dat(ch * nx + 1: (ch + 1) * nx, 1) = x;
    info_dat(ch * nx + 1: (ch + 1) * nx, 2) = t;
    info_dat(ch * nx + 1: (ch + 1) * nx, 3) = u_n;
    info_dat_analyt(ch * nx + 1: (ch + 1) * nx, 1) = x;
    info_dat_analyt(ch * nx + 1: (ch + 1) * nx, 2) = t;
    info_dat_analyt(ch * nx + 1: (ch + 1) * nx, 3) = analyt_solve;
    
    u_num_sol(ch+1, :) = u_n;
    u_analyt_sol(ch+1,:) = analyt_solve;
    
    % Solution
    t = t + dt;    
    u_nm1 = u_n;
    u_n = u_np1;
    
    % Source
%     u_n(nx / 2) = dt^2 * 600000 * cos(20*pi*t);
%     u_n(nx / 2) = dt^2 * 1000000 *  (2 * pi^2 * 40^2 * t^2) * exp(-pi^2*t^2*40^2);
%     u_n(3 * nx / 4) = dt^2 * 20000 *  (2 * pi^2 * 5^2 * t^2) * exp(-pi^2*t^2*5^2);

    for i = 2:nx - 1
        u_np1(i) = 2 * u_n(i) - u_nm1(i) +...
            CFL^2 * (u_n(i + 1) - 2 * u_n(i) + u_n(i - 1));
    end

    % Visualize
    figure(1); 
    filename = 'AnalytNumSol.gif';
    hold on;
    clf;
    plot(x, u_n, 'Color', [0.51 0.72 0.75], LW, lw);
    hold on;
    plot(x(1:1:end), analyt_solve(1:1:end), '--k', LW, lw);
    xlabel('$\xi$', 'interpreter', 'latex')
    % ylabel('$u$', 'interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment','right')
    legend('$\tilde u$', '$u$', 'interpreter', 'latex')
    title(['$\tau=$', num2str(round(t,2))], 'interpreter', 'latex')
    axis([0 Lx -0.75 0.75]);
    shg;
    pause(0.01);
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,1024);
    if ch == 0
         imwrite(imind,cm,filename,'gif', 'DelayTime',0.1, 'Loopcount',inf);
    else
         imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'WriteMode','append');
    end
    hold on;
    
    ch = ch + 1;
end

%%
% ss = '%f %f %f\n';
% FN = 'SolSinAnalyt.txt';
% fileName = fopen(FN,'w');
% fprintf(fileName, 'x t u\n');
% fprintf(fileName, ss, info_dat_analyt.');
% fclose(fileName);