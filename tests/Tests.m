% Filename: Tests.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: This script performs general tests of the geodesic utility
% functions.
% Version: 1.0

%% Clean

clear
close all
clc

%% Test 1

% Example system and metric
nx = 3;
nu = 1;

xeq = zeros(nx,1);
x = 0.01*[1;1;1];

W0 = [2.285714286 0.1428571429 -0.4285714286 ; 
      0.1428571429 1.571428571 0.2857142857;
      -0.4285714286 0.2857142857 1.142857143];
W1 = [0 -4.57142857142857 0 ; 
      -4.57142857142857 -0.571549759244019 0.857142857142854;
      0 0.857142857142854 0];
W2 = [0 0 0 ; 
      0 9.14297833067258 0;
      0 0 0];

W = @(x) (W0 + W1*x(1,:) + W2*x(1,:).^2);
dW = @(x) (W1 + 2*W2*x(1,:));

rho = @(x) (2 + 1.264*x(1,:) + 8.157*x(1,:).^2);
B = [0;0;1];

N = 4;
num_cpnts = N + 10;

% Initialize configuration (these remain always constant)
config = geodesicConfig(N,num_cpnts,3,W,dW)

% Init parameters (these remain constant during a single geodesic search)
params = geodesicParams(config, xeq, x);

opts = optimoptions('fmincon','Algorithm','sqp',...
                    'MaxFunctionEvaluations',1e6);
[copt,Eopt] = fmincon(params.f, params.c0, [], [], ...
                      params.Aeq, params.beq, [], [], [], opts);

%% Test 2 - Controller evaluation

dK = @(x) ( -0.5*rho(x)*(B'/W(x)) );

Neval = 1e3;
[u,us,epnts] = evalDiffK(config,copt,dK,false,Neval);
u

% Visualize control action along path
figure
plot(epnts,us,'LineWidth',1)
grid on
xlabel('s')
ylabel('\delta u(s)')
set(gca,'FontSize',14)
set(gcf,'Color','w')

%% Test 3 - Simulate closed-loop system

dyn = @(x,u) ([-x(1)+x(3); x(1)^2-x(2)-2*x(1)*x(3)+x(3); -x(2)+u]);

% Design LQR
Q = eye(nx);
R = eye(nu);
A = [-1, 0, 1; 0, -1, 1; 0, -1, 0];
P = are(A,B*inv(R)*B',Q);
Klqr_K = -inv(R)*B'*P;

tf = 10;               % Final time
dt = 0.05;             % Step size
x0 = 0.946*[4;4;6];          % Initial condition

% Initialize configuration (these remain always constant)
config = geodesicConfig(4, 4 + 10, 3, W, dW);
params = geodesicParams(config, xeq, x0);
xccm0 = params.c0;

[~,xs_lqr] = simRK4(dyn, {x0,[]}, tf, dt, @(x,xk) ({Klqr_K*x, []}));
[ts,xs_ccm] = simRK4(dyn, {x0,xccm0}, tf, dt, @(x,xk) ( Kccm(config,dK,xeq,x,xk) ));

figure
subplot(3,1,1)
plot(ts,xs_lqr(1,:),'LineWidth',1)
hold on
plot(ts,xs_ccm(1,:),'-r','LineWidth',1)
grid on
legend('LQR', 'CCM')
subplot(3,1,2)
plot(ts,xs_lqr(2,:),'LineWidth',1)
hold on
plot(ts,xs_ccm(2,:),'-r','LineWidth',1)
grid on
legend('LQR', 'CCM')
subplot(3,1,3)
plot(ts,xs_lqr(3,:),'LineWidth',1)
hold on
plot(ts,xs_ccm(3,:),'-r','LineWidth',1)
grid on
legend('LQR', 'CCM')
xlabel('t [s]')
