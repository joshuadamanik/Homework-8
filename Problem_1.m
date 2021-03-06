%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEWORK #8
% Joshua Julian Damanik (20194701)
% AE551 - Introduction to Optimal Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all;
addpath('lib');

eps = 0.1;
k = 15;

rho_list = [10, 100, 1000];
delta_rho = 2;

iter = 0;

X_init = [0; 10];
v_init = 0;

f = @(X) (0.5*(X(1)-1).^2+10*(X(2)-1).^2);
c = @(X) (X(1)-2)^2+2-X(2);

X = X_init;
v = v_init;

X_data = X;
v_data = v;
c_data = c(X);

for ii = 2:1000
    iter = iter + 1;
    
    L = @(X, v) (f(X) + v'*c(X));
    L_X = @(X) L(X, v);
    L_v = @(v) L(X, v);
    
    grad_L_X = grad_central_diff(X, eps, L_X);
    hess_L_X = hess_central_diff(X, eps, L_X);
    grad_c = grad_central_diff(X, eps, c);
    
    S = [hess_L_X, grad_c; grad_c', 0];
    r = [grad_L_X; c(X)];
    del_r = -pinv(S)*r;
    
    X = X + del_r(1:length(X));
    v = v + del_r(length(X)+1:end);
    
    X_data(:,ii) = X;
    v_data(:,ii) = v;
    c_data(:,ii) = c(X);
    
    fprintf('ii=%d\n',ii);
    
    if norm(X_data(:,ii)-X_data(:,ii-1))/norm(X_data(:,ii)) < 1e-10
        break;
    end
end

%% Function Graph

N_grid = 50;
axis_xy = [-2 6 0 15];
x_cont = linspace(axis_xy(1), axis_xy(2), N_grid);
y_cont = linspace(axis_xy(3), axis_xy(4), N_grid);
[X_cont, Y_cont] = meshgrid(x_cont, y_cont);

F_cont = zeros(N_grid);
Cy_cont = zeros(N_grid);

for i=1:N_grid
    for j=1:N_grid
        F_cont(i,j) = f([X_cont(i,j), Y_cont(i,j)]');
    end
end

figure(1);
s = contour(X_cont, Y_cont, F_cont);
colorbar;
hold on;

%% Constraint Graph

C_data = (x_cont-2).^2+2;
plot(x_cont, C_data, 'r--');
axis(axis_xy);

%% Search Path Graph
color = [1, 0.3, 0.3;
         0.3, 0.5, 0.3;
         0.3, 0.3, 1];
% for j=1:3
    points = X_data;
    for i=1:size(points,2)-1
        F_data = f(points(:,i));
        qlen = [points(:,i+1) - points(:,i)];% f(X_data(:,i+1))-F_data];
        quiver(points(1,i), points(2,i), ...% F_data, ...
                    qlen(1), qlen(2), ... % qlen(3), ...
                    'r', 'AutoScale', 'off', 'LineWidth', 1, ...
                    'MaxHeadSize', min(1 / norm(qlen),1), ...
                    'color', color(1,:));
    end
% end
xlabel('x');
ylabel('f(x)');
zlabel('z');
legend('Function', 'Constraint', 'Search path', 'Location', 'SouthEast');

%% v Graph
figure(2); subplot(2,1,1); hold on;
% for j=1:3
    p=plot(1:length(v_data), v_data, 'Color', color(1,:));
% end
grid on;
ylabel('v');

%% c Graph
subplot(2,1,2); hold on;
% for j=1:3
    p=plot(1:length(c_data), c_data, 'Color', color(1,:));
% end
grid on;
xlabel('Iteration');
ylabel('c');