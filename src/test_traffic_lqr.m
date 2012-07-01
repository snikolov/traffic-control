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
global A B K n_cars n_active xlast vd dd topology radius dmin tau delay

status=0;

% Gains for unactuated cars.
k1=10;
k2=10;

% Reaction Delay
delay=0;
tau=0.45;

n_cars=5;

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
  A(n_cars+1:2*n_cars,n_cars+1:2*n_cars)=k2*C1;
else
  A=zeros(n_cars*3);
  A(1:n_cars,n_cars+1:2*n_cars)=C1;
  A(n_cars+1:2*n_cars,2*n_cars+1:3*n_cars)=eye(n_cars);
  A(2*n_cars+1:3*n_cars,1:n_cars)=k1*eye(n_cars)/tau;
  A(2*n_cars+1:3*n_cars,n_cars+1:2*n_cars)=k2*C1/tau;
  A(2*n_cars+1:3*n_cars,2*n_cars+1:3*n_cars)=-eye(n_cars)/tau;
end
  
[V,L]=eig(A);

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
plot(real(diag(L)))

car_indices=1:n_cars;
active=[];
%active=car_indices(rand(n_cars,1)<0.15)

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
% Desired velocity
vd=2;
% Desired intercar distance
dd=2;
% Minimum intercar distance
dmin=0.25;

if ~delay
  x=[-1*rand(n_cars,1);vd;zeros(n_cars-1,1)];
else
  x=[-1*rand(n_cars,1);vd;zeros(n_cars-1,1);zeros(n_cars,1)];
end

if strcmpi(topology,'loop')
  radius=dd*n_cars*1.5/(2*pi);
  x(1)=radius*2*pi-sum(x(2:n_cars)+dd)-dd;
  if x(1)<-dd
    fprintf('Radius is too small');
    status=-1;
  end
end

function run(x)
global xlast vd n_cars dd dmin xmax radius tidx
% Boundary for plotting
xmax=max(x(1:n_cars));
dt=0.0001;
figure
NO_COLLIDE=0;
NO_BACKWARD=0;
for tidx=1:250000
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
  if ~mod(tidx,8000)
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
global n_cars dd vd xlast xmax topology radius tidx active

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
  title(sprintf('Cars %d',tidx))
  colormap('cool')
elseif strcmpi(topology,'loop')
  spacings=x(1:n_cars)+dd;
  thetas=pos/radius;
  theta_goals=thetas-dd/radius;
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
