% LQR-based traffic smoothing.

function test_traffic_lqr
global A B K n_cars n_active

% Gains for unactuated cars.
k1=1;
k2=1;
n_cars=5;
A=zeros(n_cars*2);
C=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
%C(1,n_cars)=1;
A(1:n_cars,n_cars+1:2*n_cars)=C;
A(n_cars+1:2*n_cars,n_cars+1:2*n_cars)=k2*C;
A(n_cars+1:2*n_cars,1:n_cars)=k1*diag(ones(n_cars,1));

% Indices of actuated cars.
active=2;
n_active=numel(active);
A(active+n_cars,:)=0;
B=eye(2*n_cars);
B=B(:,n_cars+active);

A
B

% LQR
Q=eye(2*n_cars);
R=eye(n_active);
Co=ctrb(A,B);
rank(Co)
[K,S,e]=lqr(A,B,Q,R)

% function run(x0)
% dt=0.1;
% for tidx=1:1000
%   t=dt*tidx;
%   u=control(x);
%   xdot=dynamics(x,u);
%   x=x+dt*xdot;
%   plot_cars(x)
% end
% 
% function xdot=dynamics(x,u)
% 
% function u=control(x)
% 
% function plot_cars(x)

