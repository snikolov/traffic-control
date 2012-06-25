% LQR-based traffic smoothing.

function test_traffic_lqr
global A B K n_cars n_active

% Gains for unactuated cars.
k1=1;
k2=1;
n_cars=10;
A=zeros(n_cars*2);
C=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
C(1,n_cars)=1;
A(1:n_cars,n_cars+1:2*n_cars)=C;
A(n_cars+1:2*n_cars,n_cars+1:2*n_cars)=k2*C;
A(n_cars+1:2*n_cars,1:n_cars)=k1*diag(ones(n_cars,1));

% Indices of actuated cars.
active=[1,2,3];
n_active=numel(active);
A(active+n_cars,:)=0;
B=eye(2*n_cars);
B=B(:,n_cars+active);

% LQR
Q=eye(2*n_cars);
R=eye(n_active);
Co=ctrb(A,B);
rank(Co)
[K,S,e]=lqr(A,B,Q,R)

x=[1*ones(n_cars,1);1;zeros(n_cars-1,1)];
run(x);

function run(x)
dt=0.1;
figure
for tidx=1:1000
  t=dt*tidx;
  u=control(x);
  xdot=dynamics(x,u);
  x=x+dt*xdot;
  plot_cars(x);
end

function xdot=dynamics(x,u)
global A B
xdot=A*x+B*u;

function u=control(x)
global K
u=-K*x;

function plot_cars(x)
global n_cars
pos=cumsum(x(1:n_cars));
subplot(211)
scatter(pos,zeros(size(pos)),'ks','SizeData',5)
xlim([-2*n_cars n_cars*2])
subplot(212)
stem(x(1:n_cars))
drawnow;
pause(0.01)