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
global A B K n_cars n_active dmin tau delay L topology

status=0;
topology='loop';

% Gains for unactuated cars.
k1=100;
k2=100;

n_cars=24;
L=100;
%pos=sort(rand(n_cars-1,1)*L,'descend');
dist=-rand(n_cars-1,1)*1+ones(n_cars-1,1)*L/n_cars;
x=[dist;L-sum(dist);zeros(n_cars,1)];

function run(x)
global pos1 n_cars dmin max_pos tidx L
% Boundary for plotting
max_pos=max(x(1:n_cars));
dt=0.001;
figure
pos1=0;
for tidx=1:250000
  xdot=dynamics(x);
  x=x+dt*xdot;
  x(2*n_cars)=-sum(x(n_cars+1:2*n_cars-1));
  x(n_cars)=L-sum(x(1:n_cars-1));
  pos1=pos1+dt*xdot(1);
  if ~mod(tidx,40)
    plot_cars(x,xdot);
    % disp(sum(collisions)/n_cars)
  end
end

function xdot=dynamics(x)
global n_cars
xdot=zeros(2*n_cars,1);
xdot(1:n_cars)=x(n_cars+1:2*n_cars);
xdot(n_cars+1:2*n_cars)=0.1*(vopt(x(1:n_cars))-x(n_cars+1:2*n_cars));
xdot(n_cars)=-sum(xdot(1:n_cars-1));
xdot(2*n_cars)=-sum(xdot(n_cars+1:2*n_cars-1));

function v=vopt(h)
v=zeros(size(h));
v(h>1)=5*(h(h>1)-1).^3./(1+(h(h>1)-1)).^3;

function plot_cars(x,xdot)
global n_cars pos1 max_pos topology radius tidx active

pos=[pos1;pos1-cumsum(x(2:n_cars))];

subplot(4,2,[1,3,5,7])
if strcmpi(topology,'line')
  if pos1>max_pos
    max_pos=pos1+15;
  end
  scatter(pos,zeros(size(pos)),linspace(10,80,n_cars),linspace(1,32,n_cars),'filled')
  xlim([pos1-100*n_cars*1.0, max_pos])
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
