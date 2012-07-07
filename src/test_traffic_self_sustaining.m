% LQR-based traffic smoothing.

function test_traffic_lqr
global a
%rng('default')

close all

count=0;
%as=linspace(1,2,3);
while count<1%numel(as)
  %a=as(count+1);
  [x,status]=init;
  if status~=-1
    run(x)
    count=count+1;
  end
end

% init;
% roa;

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function [x,status]=init
global n_cars L a vmax active
status=0;
a=1.5;

% rng('default');

vmax=1;
n_cars=200;
L=400;

pos=flipud(linspace(L/n_cars,L,n_cars)');
pert=0.1*L/n_cars*(2*rand(n_cars,1)-1);
%pert=0.3*L/n_cars*(rand(n_cars,1)<0.01);
pos=pos+pert;

vel=vopt(xtod(pos));
%pert=0.3*L/n_cars*(rand(n_cars,1)<0.01);
vel=vel+pert;
x=[pos;vel];
%x=[pos;zeros(n_cars,1)];
%x(floor(1.1*n_cars))=1;

%car_ind=2:n_cars;
%active=car_ind(rand(1,n_cars-1)<0.1);
%active=[1:50:n_cars];
active=[];

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function run(x)
global n_cars tidx L a
% Boundary for plotting
dt=0.3;
figure
iter=2000;
skip=10;
time_evol=zeros(n_cars*floor(iter/skip),2);
for tidx=1:iter
  xdot=dynamics(x);
  if ~mod(tidx,skip)
    %plot_cars(x,xdot);
    time_evol(n_cars*(tidx/skip-1)+1:n_cars*(tidx/skip),:)=horzcat(mod(x(1:n_cars),L),tidx*dt*ones(n_cars,1));
  end
  
  % Runge-Kutta integration.
  xdotk1=dynamics(x);
  xt=x+dt*xdotk1*0.5;
  
  xdotk2=dynamics(xt);
  xt=x+xdotk2*dt*0.5;
  
  xdotk3=dynamics(xt);
  xt=x+xdotk3*dt;
  
  xdotk4=dynamics(xt);
  x=x+(xdotk1+2*xdotk2+2*xdotk3+xdotk4)/6*dt;
end

scatter(time_evol(:,1),time_evol(:,2),'k.','SizeData',1);
title(sprintf('%.4f',a));
set(gcf,'Position',[200,200,400,300]);

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function xdot=dynamics(x)
global n_cars a active
xdot=zeros(2*n_cars,1);
d=xtod(x);
xdot(1:n_cars)=x(n_cars+1:2*n_cars);
xdot(n_cars+1:2*n_cars)=a*(vopt(d)-x(n_cars+1:2*n_cars));
% Linear PD
% xdot(active+n_cars)=0.001*(x(active-1)-x(active)-1*x(n_cars+active))+0.001*(x(active+n_cars-1)-x(active+n_cars));
% Nonlinear: Pretend there is less headway (be more conservative)
xdot(n_cars+active)=a*(vopt(d(active).^0.5)-x(n_cars+active));

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function d=xtod(x)
global L n_cars
d=[x(n_cars)-x(1)+L;x(1:n_cars-1)-x(2:n_cars)];

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function v=vopt(h)
global vmax L n_cars
v=zeros(size(h));
%v(h>1)=vmax*(h(h>1)-1).^3./(1+(h(h>1)-1).^3);
%v(h>1)=vmax*(h(h>1)-1).^3./h(h>1).^3;
v=vmax*(tanh(h-2)+tanh(2));

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function roa
global n_cars a
f=1;
A=zeros(n_cars*2);
A(1,1+n_cars)=1;
A(1+n_cars,n_cars)=a*f;
A(1+n_cars,1)=-a*f;
A(1+n_cars,1+n_cars)=-a;
for i=2:n_cars
  A(i,i+n_cars)=1;
  A(i+n_cars,i-1)=a*f;
  A(i+n_cars,i)=-a*f;
  A(i+n_cars,i+n_cars)=-a;
end

P=lyap(A,eye(2*n_cars));
max(P(:))
res=50;
[X,Y]=meshgrid(linspace(-13,13,res),linspace(-13,13,res));
x=[reshape(X,1,numel(X));reshape(Y,1,numel(Y))];
for i=1:50;
  l=ceil(rand*2*n_cars);
  m=ceil(rand*2*n_cars);
  V=diag(x'*P([l,m],[l,m])*x);
  V=reshape(V,res,res);
  imagesc(V)
  colorbar;
  pause;
end


%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function plot_cars(x,xdot)
global n_cars L tidx active

subplot(4,2,[1,3,5,7])
  radius=L/(2*pi);
  thetas=x(1:n_cars)/radius;
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

subplot(422)
%stem(xtod(x(1:n_cars)),'k')
plot(xtod(x(1:n_cars)),'k')
title('Spacings')

subplot(424)
%stem(x(1:n_cars),'k')
plot(x(1:n_cars),'k')
title('Positions')

subplot(426)
%stem(x(n_cars+1:2*n_cars),'k')
plot(x(n_cars+1:2*n_cars),'k')
title('Velocities')

subplot(428)
%stem(xdot(n_cars+1:2*n_cars),'k')
plot(xdot(n_cars+1:2*n_cars),'k')
title('Accelerations')

drawnow;
set(gcf,'Position',[100,100,1000,500]);
pause(0.01)
