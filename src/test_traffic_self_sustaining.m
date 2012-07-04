% LQR-based traffic smoothing.

function test_traffic_lqr

rng('default')

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
global n_cars L a vmax
status=0;

a=5;
vmax=5;
n_cars=15;
L=50;
pos=sort(rand(n_cars,1)*L,'descend');
x=[pos;zeros(n_cars,1)];

function run(x)
global pos n_cars max_pos tidx L
% Boundary for plotting
max_pos=max(x(1:n_cars));
dt=0.05;
figure
pos=[0;-cumsum(x(2:n_cars))];
for tidx=1:250000
  xdot=dynamics(x);
  x=x+dt*xdot;
  if ~mod(tidx,5)
    plot_cars(x,xdot);
    % disp(sum(collisions)/n_cars)
  end
end

function xdot=dynamics(x)
global n_cars a
xdot=zeros(2*n_cars,1);
d=xtod(x);
xdot(1:n_cars)=x(n_cars+1:2*n_cars);
xdot(n_cars+1:2*n_cars)=a*(vopt(d)-x(n_cars+1:2*n_cars));

function d=xtod(x)
global L n_cars
d=[x(n_cars)-x(1)+L;x(1:n_cars-1)-x(2:n_cars)];

function v=vopt(h)
global vmax
v=zeros(size(h));
v(h>1)=vmax*(h(h>1)-1).^3./(1+(h(h>1)-1)).^3;
%v(h>1)=vmax*(tanh(h(h>1)-2)+tanh(2));

function plot_cars(x,xdot)
global n_cars L tidx active

subplot(4,2,[1,3,5,7])
  radius=L/(2*pi);
  thetas=x(1:n_cars)/radius;
  scatter(radius*cos(thetas),radius*sin(thetas),linspace(10,80,n_cars),'r','filled');
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

subplot(422)
stem(xtod(x(1:n_cars)))
title('Spacings')

subplot(424)
stem(x(1:n_cars))
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
