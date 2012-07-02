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
global A B K n_cars active xlast vd dd topology radius dmin tau delay

status=0;

% Gains for unactuated cars.
k1=15;
k2=15;

% Reaction Delay
delay=0;
tau=0.1;

n_cars=25;

% Build the system matrices.
topology='line';
C1=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
C2=diag(-1*ones(n_cars,1));%+diag(ones(n_cars-1,1),-1);
if strcmpi(topology,'loop')
  C1(1,n_cars)=1;
end

if ~delay
  A=zeros(n_cars*2);
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
%active=car_indices(rand(n_cars,1)<0.25);

n_active=numel(active);
if n_active>0
  if ~delay
    A(active+n_cars,:)=0;
    B=eye(2*n_cars);
    B=B(:,n_cars+active);
    Q=eye(2*n_cars);
    Q(1:n_cars+1,1:n_cars+1)=eye(n_cars+1);
  else
    A(active+2*n_cars,1:2*n_cars)=0;
    B=eye(3*n_cars)/tau;
    B=B(:,2*n_cars+active);
    Q=zeros(3*n_cars);
    Q(1:n_cars+1,1:n_cars+1)=eye(n_cars+1);
  end
  R=eye(n_active);
  try
    [K,S,e]=lqr(A,B,Q,R)
  catch exc
    fprintf('LQR did not work\n');
    status=-1;
  end
else
  B=0;
  K=0;
end

% Keep track of the position of the last car
xlast=0;
% Desired velocity
vd=1;
% Desired intercar distance
dd=5;
% Minimum intercar distance
dmin=2;

if ~delay
  x=[dd*rand(n_cars,1)-dd/2;zeros(n_cars,1)];
else
  x=[dd*rand(n_cars,1)-dd/2;zeros(n_cars,1);zeros(n_cars,1)];
end

if strcmpi(topology,'loop')
  radius=dd*n_cars*0.95/(2*pi);
  x(1)=radius*2*pi-sum(x(2:n_cars)+dd)-dd;
  if x(1)<-dd
    fprintf('Radius is too small');
    status=-1;
  end
end

function run(x)
global xlast vd n_cars dd dmin xmax tidx delay
% Boundary for plotting
xmax=max(x(1:n_cars));
dt=0.01;
figure
NO_COLLIDE=0;
NO_BACKWARD=0;
for tidx=1:250000
  u=control(x);
  xdot=dynamics(x,u);
  if NO_BACKWARD
    xdot(1:n_cars)=max(xdot(1:n_cars),-vd);
    x(n_cars+1:2*n_cars)=max(x(n_cars+1:2*n_cars),-vd);
  end
  %x(1)=2*pi*radius-sum(x(2:n_cars)+dd)-dd;
  if NO_COLLIDE
    for cidx=1:n_cars
      x(cidx)=x(cidx)+xdot(cidx)*dt;
      if x(cidx)<dmin-dd;
        x(cidx+n_cars)=-vd;
        x(cidx)=dmin-dd;
      end
    end
    x(n_cars+1:2*n_cars)=x(n_cars+1:2*n_cars)+dt*xdot(n_cars+1:2*n_cars); 
    if delay
      x(2*n_cars+1:end)=x(2*n_cars+1:end)+dt*xdot(2*n_cars+1:end);    
    end
%     delta=0;
%     collision=1;
%     while collision
%       collision=0;
%       for cidx=1:n_cars
%         x(cidx)=x(cidx)-delta;
%         if x(cidx)<dmin-dd
%           x(cidx)=dmin-dd;
%           x(cidx+n_cars)=0;
%           delta=dmin-dd-x(cidx)
%           collision=1;
%           %plot_cars(x,xdot);
%           %pause;
%         else
%           delta=0;
%         end
%       end
%     end
    
  else
    x=x+dt*xdot;
  end
  xlast=xlast+dt*(vd+xdot(2*n_cars));
  if ~mod(tidx,100)
    plot_cars(x,xdot);
    % disp(sum(collisions)/n_cars)
  end
end

function xdot=dynamics(x,u)
global A B
xdot=A*x+B*u;

function u=control(x)
global K active
if active
  xhat=estimate_x(x);
  u=-K*xhat;
else
  u=-K*x;
end
  
function xhat=estimate_x(x)
global active delay n_cars
xavg=mean(x(active));
vavg=mean(x(n_cars+active));
xhat=[xavg*ones(n_cars,1);vavg*ones(n_cars,1)];
xhat(active)=x(active);
xhat(active+n_cars)=x(active+n_cars);
if delay
  aavg=mean(x(2*n_cars+active));
  xhat=[xhat;aavg*ones(n_cars,1)];
  xhat(active+2*n_cars)=x(active+2*n_cars);
end


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
