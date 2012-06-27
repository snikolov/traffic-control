% LQR-based traffic smoothing.

function test_traffic_lqr
global A B K n_cars n_active x1 vd dd topology

close all

% Gains for unactuated cars.
k1=2.15;
k2=2.15;
n_cars=30;
A=zeros(n_cars*2);
C=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
topology='loop';
if strcmpi(topology,'loop')
  C(1,n_cars)=1;
end
A(1:n_cars,n_cars+1:2*n_cars)=C;
A(n_cars+1:2*n_cars,n_cars+1:2*n_cars)=k2*C;
A(n_cars+1:2*n_cars,1:n_cars)=k1*diag(ones(n_cars,1));
if strcmpi(topology,'line')
  A(1,n_cars+1)=0;
  A(n_cars+1,1)=0;
end

% Indices of actuated cars.
active=[]
n_active=numel(active);
if n_active>0
  A(active+n_cars,:)=0;
  B=eye(2*n_cars);
  B=B(:,n_cars+active);
  % LQR
  Q=eye(2*n_cars);
  R=eye(n_active);
  Co=ctrb(A,B);
  rank(Co)
  [K,S,e]=lqr(A,B,Q,R)
else
  B=0;
  K=0;
end

% Keep track of the position of the first car
x1=0;
% Desired velocity
vd=10;
% Desired intercar distance
dd=5;

x=[0.1*rand(n_cars,1);vd;zeros(n_cars-1,1)];
run(x);

function run(x)
global x1 vd n_cars dd xmax
% Boundary for plotting
xmax=max(x(1:n_cars));
dt=0.1;
figure
NO_COLLIDE=0;
NO_BACKWARD=0;
for tidx=1:10000
  u=control(x);
  xdot=dynamics(x,u);
  if ~NO_COLLIDE
    % Make sure cars don't collide.
    collisions=x(1:n_cars)<-dd-dd/10;
    x(collisions)=0;
    xdot(collisions)=0;
  end
  if ~NO_BACKWARD
    % Make sure cars don't move backward.
    backward=x(n_cars+1:2*n_cars)<0;
    x(n_cars+backward)=0;
    xdot(backward)=0;
  end
  x=x+dt*xdot;
  x1=x1+dt*(vd+xdot(n_cars+1));
  if ~mod(tidx,5)
    plot_cars(x,xdot);
    if x(1)<-dd/2
      1
    end
    % disp(sum(collisions)/n_cars)
  end
end

function xdot=dynamics(x,u)
global A B
xdot=A*x+B*u;

function u=control(x)
global K
u=-K*x;

function plot_cars(x,xdot)
global n_cars dd vd x1 xmax topology

pos=[x1;x1-(cumsum(x(2:n_cars)+dd))];

subplot(4,2,[1,3,5,7])
if strcmpi(topology,'line')
  if x1>xmax
    xmax=x1+5*dd;
  end
  scatter(pos,zeros(size(pos)),'ks','SizeData',20)
  xlim([x1-1.6*n_cars*dd, xmax])
  ylim([-2,2])
  title('Cars')
elseif strcmpi(topology,'loop')
  spacings=x(1:n_cars)+dd;
  radius=sum(spacings(2:end))/(2*pi);
  theta1=x1/radius;
  thetas=[theta1;theta1-cumsum(spacings(2:end))/radius];
  scatter(radius*cos(thetas),radius*sin(thetas),linspace(10,80,n_cars),linspace(1,32,n_cars),'filled');
  colormap('cool')
  axis square
  xlim([-radius,radius])
  ylim([-radius,radius])
end

subplot(422)
stem(dd+x(1:n_cars))
title('Spacings')

subplot(424)
stem(pos)
title('Positions')

subplot(426)
stem(x(vd+n_cars+1:2*n_cars))
title('Velocities')

subplot(428)
stem(xdot(n_cars+1:2*n_cars))
title('Accelerations')

drawnow;
set(gcf,'Position',[100,100,1000,500]);
pause(0.01)
