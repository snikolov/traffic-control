% LQR-based traffic smoothing.

function test_traffic_lqr

close all
count=0;
while count<10
  [x,status]=init;
  if status~=-1
    run(x)
    count=count+1;
  end
end

function [x,status]=init
global A B K n_cars n_active xlast topology radius dmin tau delay

status=0;

% Gains for unactuated cars.
k1=100;
k2=100;

% Reaction Delay
delay=1;
tau=0.35;

n_cars=125;
L=10;

% Build the system matrices.
C1=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
C2=diag(-1*ones(n_cars,1));%+diag(ones(n_cars-1,1),-1);
if strcmpi(topology,'loop')
  C1(1,n_cars)=1;
end

if ~delay
  A=zeros(n_cars*2);
  topology='loop';
  A(1:n_cars,n_cars+1:2*n_cars)=C1;
  A(n_cars+1:2*n_cars,1:n_cars)=k1*eye(n_cars);
  A(n_cars+1:2*n_cars,n_cars+1:2*n_cars)=k1*C2;
else
  A=zeros(n_cars*3);
  A(1:n_cars,n_cars+1:2*n_cars)=C1;
  A(n_cars+1:2*n_cars,2*n_cars+1:3*n_cars)=eye(n_cars);
  A(2*n_cars+1:3*n_cars,1:n_cars)=k1*eye(n_cars)/tau;
  A(2*n_cars+1:3*n_cars,n_cars+1:2*n_cars)=k1*C2/tau;
  A(2*n_cars+1:3*n_cars,2*n_cars+1:3*n_cars)=-eye(n_cars)/tau;
end
  
[V,Lm]=eig(A);

%if strcmpi(topology,'line')
%  A(1,n_cars+1)=0;
%  A(n_cars+1,1)=0;
%end

% Indices of actuated cars.

figure(45)
subplot(211)
imagesc(abs(V))
colorbar
subplot(212)
plot(real(diag(Lm)))

car_indices=1:n_cars;
active=[];
%active=car_indices(rand(n_cars,1)<0.15)

n_active=numel(active);
if n_active>0
  if ~delay
    A(active+n_cars,:)=0;
    B=eye(2*n_cars);
    B=B(:,n_cars+active);
  else
    A(active+2*n_cars,1:2*n_cars)=0;
    B=eye(3*n_cars)/tau;
    B=B(:,2*n_cars+active);
  end
  
  % LQR
  Q=eye(2*n_cars);
  R=eye(n_active);
  Co=ctrb(A,B);
  rank(Co)
  try
    [K,S,e]=lqr(A,B,Q,R)
  catch exc
    fprintf('LQR did not work');
    status=-1;
  end
else
  B=0;
  K=0;
end

% Keep track of the position of the last car
xlast=0;
% Minimum intercar distance
dmin=0.25;

if ~delay
  x=[L*rand(n_cars,1)/n_cars;1;zeros(n_cars-1,1)];
else
  x=[L*rand(n_cars,1)/n_cars;1;zeros(n_cars-1,1);zeros(n_cars,1)];
end

if strcmpi(topology,'loop')
  radius=L/(2*pi);
  x(1)=radius*2*pi-sum(x(2:n_cars));
end

function run(x)
global xlast n_cars dmin xmax radius tidx
% Boundary for plotting
xmax=max(x(1:n_cars));
dt=0.001;
figure
NO_COLLIDE=0;
NO_BACKWARD=0;
for tidx=1:250000
  u=control(x);
  xdot=dynamics(x,u);
  x=x+dt*xdot;
  x(1)=2*pi*radius-sum(x(2:n_cars));
  if NO_COLLIDE
    % Make sure cars don't collide.
    collisions=1+x(2:n_cars)<dmin;
    x(collisions)=dmin;
    x(n_cars+collisions)=0;
    xdot(collisions)=0;
  end
  if NO_BACKWARD
    xdot(1:n_cars)=max(xdot(1:n_cars),0);
  end
  xlast=xlast+dt*(xdot(2*n_cars));
  if ~mod(tidx,4000)
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
global n_cars xlast xmax topology radius tidx active

pos=flipud([xlast;xlast+cumsum(flipud(x(2:n_cars)))]);

subplot(4,2,[1,3,5,7])
if strcmpi(topology,'line')
  x1=xlast+sum(x(2:n_cars));
  if x1>xmax
    xmax=x1+15;
  end
  scatter(pos,zeros(size(pos)),linspace(10,80,n_cars),linspace(1,32,n_cars),'filled')
  xlim([x1-100*n_cars*1.0, xmax])
  ylim([-2,2])
  title(sprintf('Cars %d',tidx))
  colormap('cool')
elseif strcmpi(topology,'loop')
  spacings=x(1:n_cars);
  thetas=pos/radius;
  scatter(radius*cos(thetas),radius*sin(thetas),linspace(10,80,n_cars),linspace(1,32,n_cars),'filled');
  colormap('cool')
  hold on
  for aidx=active
    scatter(radius*cos(thetas(aidx)),radius*sin(thetas(aidx)),'ro','SizeData',25,'MarkerEdgeColor','r','MarkerFaceColor','r');
  end
  hold off
  % hold on
  % scatter(radius*cos(theta_goals),radius*sin(theta_goals),5,'k','filled');
  hold off
  axis square
  xlim([-radius,radius])
  ylim([-radius,radius])
  title(sprintf('Cars %d',tidx))
end

subplot(422)
stem(x(1:n_cars))
title('Spacings')

subplot(424)
stem(pos)
title('Positions')

subplot(426)
stem(x(n_cars+1:2*n_cars))
title('Velocities')

subplot(428)
stem(xdot(n_cars+1:2*n_cars))
title('Accelerations')

drawnow;
set(gcf,'Position',[100,100,1000,500]);
pause(0.01)
