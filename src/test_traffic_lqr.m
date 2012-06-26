% LQR-based traffic smoothing.

function test_traffic_lqr
global A B K n_cars n_active x1 vd dd

% Gains for unactuated cars.
k1=2;
k2=2;
n_cars=55;
A=zeros(n_cars*2);
C=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
%C(1,n_cars)=1;
A(1:n_cars,n_cars+1:2*n_cars)=C;
A(n_cars+1:2*n_cars,n_cars+1:2*n_cars)=k2*C;
A(n_cars+1:2*n_cars,1:n_cars)=k1*diag(ones(n_cars,1));

% Indices of actuated cars.
active=[14,15]
n_active=numel(active);
A(active+n_cars,:)=0;
B=eye(2*n_cars);
B=B(:,n_cars+active);

% Keep track of the position of the first car
x1=0;
% Desired velocity
vd=2;
% Desired intercar distance
dd=2;

% LQR
Q=eye(2*n_cars);
R=eye(n_active);
Co=ctrb(A,B);
rank(Co)
[K,S,e]=lqr(A,B,Q,R)

x=[1*rand(n_cars,1);1;zeros(n_cars-1,1)];
run(x);

function run(x)
global x1 vd n_cars dd xmax
% Boundary for plotting
xmax=max(x(1:n_cars));
dt=0.01;
figure
for tidx=1:100000
  u=control(x);
  xdot=dynamics(x,u);
  collisions=x(1:n_cars)<-dd/10;
  x(collisions)=0;
  xdot(collisions)=0;
  x=x+dt*xdot;
  x1=x1+dt*(vd+xdot(n_cars+1));
  if ~mod(tidx,100)
    plot_cars(x);
    % disp(sum(collisions)/n_cars)
  end
end

function xdot=dynamics(x,u)
global A B
xdot=A*x+B*u;

function u=control(x)
global K
u=-K*x;

function plot_cars(x)
global n_cars dd vd x1 xmax 

pos=x1-(cumsum(x(1:n_cars)+dd));

subplot(411)
scatter(pos,zeros(size(pos)),'ks','SizeData',10)
if x1>xmax
  xmax=x1+5*dd;
end
xlim([x1-dd*n_cars*dd, xmax])
ylim([-2,2])
title('Cars')

subplot(412)
stem(x(1:n_cars)+dd)
title('Spacings')

subplot(413)
stem(pos)
title('Positions')

subplot(414)
stem(vd+x(n_cars+1:2*n_cars))
title('Velocities')
drawnow;
pause(0.01)
