%{
 =======================================================
| Program to implement three explicit finite-difference |
| methods to solve 1D linear transport equation         |
 -------------------------------------------------------
| A. 1st Order Upwind Method                            |
| B. Lax-Wendroff Method                                |
| C. 4th order Runge-Kutta                              |
 -------------------------------------------------------
| The solutions from the numerical methods can then be  |
| compared with analytical solutions.                   |
 =======================================================

 =======================================================
| Variables                                             |
 -----------                                            |
| x_min:  lower boundary of the spatial domain, meters  |
| x_max:  upper boundary of the spatial domain, meters  |
| x:      spatial domain, meters                        |
| t:      time, seconds                                 |
| c:      characteristic velocity, meters/second        |
| c_bar:  Courant number                                |
| mu:     mean                                          |
| sigma:  standard deviation                            |
| phi_exact:  exact solution                            |
 =======================================================
%}

clc;    clear;    close all;

% describe the spatial dimension (x-coordinate)
x_min = 0;
x_max = 100;
N = 200;
% vector x
dx = (x_max-x_min)/N;
x = x_min:dx:x_max;

%{
 =======================================================
| Evaluate solution for two different Courant numbers   |
| A) Close to the stability limit: 0.999                |
| B) Below the stability limit: 0.5                     |
 -------------------------------------------------------
| May run the program  by changing values of Courant    |
| number to study the solutions at various values of    | 
| c_bar.                                                |
 =======================================================
%}
c = .5;         % wave speed
c_bar = .5;     % Courant's number

t_max = 100;
dt = c_bar*dx/c;
t = 0:dt:t_max;

%{
 =======================================================
| Initial Condition                                     |
 -------------------                                    |
| Gaussian distribution with mean 0 and sd 1            |
 =======================================================
%}
mu = 0;
sigma = 1;

% Initial conditions
phi = zeros(length(x), length(t));
phi(1,:) = 0;
phi(:,1) = 1/(sigma * sqrt(2*pi)) * exp(-1/2*((x- mu)/sigma).^2);


% Calculate exact/analytical solution
phi_exact = 1/(sigma * sqrt(2*pi)) * exp(-1/2*((x -c*t_max - mu)/sigma).^2);
plot(x, phi_exact, "LineWidth", 1.5)
hold on

% 1st order Upwind method
phi1 = upwind(x, t, phi, c_bar);
plot(x, phi1(:,t_max), 'k-.', "LineWidth", 1.5);
hold on

% Lax-Wendroff method
phi2 = lax_wendroff(x,t,phi,c_bar);
plot(x, phi2(:, t_max), 'r--', "LineWidth", 1.5);

% 4th order Runge-Kutta (rk4)
phi3 = rk4(x, t, phi, dx, dt, c);
plot(x, phi2(:, t_max), 'g:', "LineWidth", 2.5);

% Set labels, title, and legend
xlabel("x");
ylabel("\phi")
head = sprintf("\\phi vs x at t = %.0f and Courant %.3f", t_max,c_bar);
title(head)
dbclear all
legend("Exact", "Upwind", "Lax-Wendroff", "Runge-Kutta");

% Comparison of Analytical solution and numerical solution
% on a 3-dimensional plot

figure(2)
[xx, tt] = meshgrid(x, t);

% 3-dimensional exact solution
phi_ex_3d = 1/(sigma * sqrt(2*pi)) * exp(-1/2*((xx -c*tt - mu)/sigma).^2);
plot3(xx, tt, phi_ex_3d, 'r')
xlabel("x");
ylabel("Time, t (s)");
zlabel("\phi")
hold on;

% Can select the method of choice 
% Edit title accordingly
plot3(xx, tt, phi2, 'k'); legend("Exact")
a = sprintf("Comparison of Exact Solution and Lax-Wendroff method until t = %.f, and Courant %.3f", t_max, c_bar);
title(a)
