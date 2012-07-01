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
global A B K n_cars n_active x1 topology radius dmin tau delay

status=0;

% Gains for unactuated cars.
k1=100;
k2=100;

% Reaction Delay
delay=0;
tau=0.35;

n_cars=35;
L=100;

topology='loop';
% Build the system matrices.
C1=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
C2=diag(-1*ones(n_cars,1));%+diag(ones(n_cars-1,1),-1);
if strcmpi(topology,'loop')
  C1(1,n_cars)=1;
end

if ~delay
  A=zeros(n_cars*2);
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

% Keep track of the position of the first car
x1=0;
% Minimum intercar distance
dmin=0.25;

pos=sort(rand(n_cars,1)*L,'descend');
dist=pos(1:n_cars-1)-pos(2:n_cars);
x=[L-sum(dist);dist;0;zeros(n_cars-1,1)];
if delay
  x=[x;zeros(n_cars,1)];
end

radius=L/(2*pi);

function run(x)
global x1 n_cars dmin xmax radius tidx
% Boundary for plotting
xmax=max(x(1:n_cars));
dt=0.001;
figure
NO_COLLIDE=0;
NO_BACKWARD=0;
for tidx=1:250000
  xdot=dynamics(x);
  x=x+dt*xdot;
  x(1)=2*pi*radius-sum(x(2:n_cars));
  if NO_COLLIDE
    % TODO
  end
  if NO_BACKWARD
    xdot(1:n_cars)=max(xdot(1:n_cars),0);
  end
  x1=x1+dt*xdot(1);
  if ~mod(tidx,40)
    plot_cars(x,xdot);
    % disp(sum(collisions)/n_cars)
  end
end

function xdot=dynamics(x)
global A B n_cars
xdot=zeros(2*n_cars,1);
xdot(1:n_cars)=x(n_cars+1:2*n_cars);
xdot(n_cars+1:2*n_cars)=1*(vopt(x(1:n_cars))-x(n_cars+1:2*n_cars));

function v=vopt(h)
v=zeros(size(h));
v(h>1)=1*(h(h>1)-1).^3./(1+(h(h>1)-1)).^3;

function plot_cars(x,xdot)
global n_cars x1 xmax topology radius tidx active

pos=[x1;x1-cumsum(x(2:n_cars))];

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
