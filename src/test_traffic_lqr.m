% LQR-based traffic smoothing.

function test_traffic_lqr
global A B K n_cars n_active xlast vd dd topology radius dmin

close all
rng('default')

% Gains for unactuated cars.
k1=2;
k2=2;
n_cars=25;
A=zeros(n_cars*2);
C1=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
C2=diag(-1*ones(n_cars,1));%+diag(ones(n_cars-1,1),-1);
topology='loop';
if strcmpi(topology,'loop')
  C1(1,n_cars)=1;
end
A(1:n_cars,n_cars+1:2*n_cars)=C1;
A(n_cars+1:2*n_cars,n_cars+1:2*n_cars)=k2*C1;
A(n_cars+1:2*n_cars,1:n_cars)=k1*diag(ones(n_cars,1));

%if strcmpi(topology,'line')
%  A(1,n_cars+1)=0;
%  A(n_cars+1,1)=0;
%end

% Indices of actuated cars.
active=[5,6]
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

% Keep track of the position of the last car
xlast=0;
% Desired velocity
vd=2;
% Desired intercar distance
dd=2;
% Minimum intercar distance
dmin=0.25;

x=[1*rand(n_cars,1)-0.5;vd;zeros(n_cars-1,1)];
if strcmpi(topology,'loop')
  radius=dd*n_cars*1.05/(2*pi);
  x(1)=radius*2*pi-sum(x(2:n_cars)+dd)-dd;
  if x(1)<-dd
    error('Radius is too small');
  end
end
run(x);

function run(x)
global xlast vd n_cars dd dmin xmax radius
% Boundary for plotting
xmax=max(x(1:n_cars));
dt=0.01;
figure
NO_COLLIDE=0;
NO_BACKWARD=0;
for tidx=1:100000
  u=control(x);
  xdot=dynamics(x,u);
  x=x+dt*xdot;
  if NO_COLLIDE
    % Make sure cars don't collide.
    collisions=1+x(2:n_cars)<dmin-dd;
    x(collisions)=dmin-dd;
    x(n_cars+collisions)=-vd;
    xdot(collisions)=0;
  end
  if NO_BACKWARD
%     % Make sure cars don't move backward.
%     backward=x(n_cars+1:2*n_cars)<0;
%     x(n_cars+backward)=-vd;
%     xdot(backward)=-vd;
%     xdot(n_cars+backward)=0;
    xdot(1:n_cars)=max(xdot(1:n_cars),-vd);
  end
  %x(1)=max(dmin,2*pi*radius-sum(x(2:n_cars)+dd))-dd;
  xlast=xlast+dt*(vd+xdot(2*n_cars));
  if ~mod(tidx,5)
    plot_cars(x,xdot);
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
global n_cars dd vd xlast xmax topology radius

pos=flipud([xlast;xlast+cumsum(flipud(x(2:n_cars))+dd)]);

subplot(4,2,[1,3,5,7])
if strcmpi(topology,'line')
  x1=xlast+sum(x(2:n_cars)+dd);
  if x1>xmax
    xmax=x1+5*dd;
  end
  scatter(pos,zeros(size(pos)),linspace(10,80,n_cars),linspace(1,32,n_cars),'filled')
  xlim([x1-dd*n_cars*1.0, xmax])
  ylim([-2,2])
  title('Cars')
  colormap('cool')
elseif strcmpi(topology,'loop')
  spacings=x(1:n_cars)+dd;
  thetas=pos/radius;
  theta_goals=thetas-dd/radius;
  scatter(radius*cos(thetas),radius*sin(thetas),linspace(10,80,n_cars),linspace(1,32,n_cars),'filled');
  % hold on
  % scatter(radius*cos(theta_goals),radius*sin(theta_goals),5,'k','filled');
  hold off
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
stem(vd+x(n_cars+1:2*n_cars))
title('Velocities')

subplot(428)
stem(xdot(n_cars+1:2*n_cars))
title('Accelerations')

drawnow;
set(gcf,'Position',[100,100,1000,500]);
pause(0.01)
