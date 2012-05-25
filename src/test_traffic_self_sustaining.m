%test_traffic
function test_traffic_self_sustaining
global k_pd d_f xmax active_cars dt tidx xdothist

close all

k_pd=1;
d_f=0.005;

fig_handle=345;
figure(fig_handle);
% Forward simulate dynamics.
dt=0.0001;
n_steps=10000;
n_cars=100;
active_cars=[];%4:10:70;
T=0:dt:dt*n_steps;
d_init=1*d_f;
% Positions.
% x=[linspace(n_cars*d_init,d_init,n_cars)';[1;zeros(n_cars-1,1)]];
%x=15*rand(n_cars,1);
%x=sort(x,'descend');
x=linspace(1.1+(n_cars-1)*d_init,1.1,n_cars)';
xdothist=[];
% Velocities.
x=[x;1;zeros(n_cars-1,1)];
%x=[x;0.4+0.05*rand(n_cars,1)];
xmax=max(x(1:n_cars));

for tidx=1:numel(T)
  t=T(tidx);
  xdot=dynamics(x,t,0);
%   x(n_cars+1:end)=x(n_cars+1:end)+dt*xdot(n_cars+1:end);
%   if x(1)+dt*xdot(1)>x(n_cars)+xmax-d_f
%     x(1)=x(n_cars)+xmax-d_f;
%     %x(1+n_cars)=0;
%   end
%   for i=2:n_cars
%     if x(i)+dt*xdot(i)>x(i-1)-d_f;
%       x(i)=x(i-1)-d_f;
%       %x(i+n_cars)=0;    
%     end
%   end
%   x=x+dt*xdot;

  % Throwaway x_future to report collisions before handling them.
%   x_future=x+dt*xdot;
%   if any((x_future(1:n_cars-1)-x_future(2:n_cars))<0)
%     fprintf('COLLISION!\n');
%   end
  
  x(1)=min(x(1)+dt*xdot(1),x(n_cars)+xmax-d_f);
  for i=2:n_cars
    x(i)=min(x(i)+dt*xdot(i),x(i-1)-d_f);
  end
  x(n_cars+1:end)=x(n_cars+1:end)+dt*xdot(n_cars+1:end);
end

function xdot=dynamics(x,t,u)
global xmax active_cars k_pd d_f dt xdothist
n_cars=numel(x)/2;
xbar1=(x(n_cars)+xmax)-x(1);
xbar=[xbar1;x(1:n_cars-1)-x(2:n_cars)];
xdot=zeros(2*n_cars,1);
vmax=25;
xdot(1:n_cars)=vmax*(1-d_f./xbar);

xdothist=[xdothist,xdot];
% Delay
memory=100;
tau=50;
weights=exp(-(1:1:memory)'/tau);
weights=weights/sum(weights);
if size(xdothist,2)>=memory
  xdot=xdothist*weights;
  xdothist=xdothist(:,2:end);
end


% % Dumb cars
% T=50*dt;
% c=0.5;
% dumb_cars=logical(zeros(n_cars,1));
% dumb_cars(10)=1;
% %modulator=double(mod(t,T)>c*T);
% modulator=(1+sin(2*pi*t/T))/2;
% xdot(dumb_cars)=xdot(dumb_cars)*modulator;

%Active cars
for i=active_cars
%   % Bando heuristic.
%   front=active_cars(active_cars>=i);
%   back=active_cars(active_cars<i);
%   xdot_bando=0;
%   bando_weights=0;
%   for j=front
%     dist=x(j)-x(i);
%     weight=exp(-dist/(n_cars*d_f/5));
%     bando_weights=bando_weights+weight;
%     xdot_bando=xdot_bando-10*weight*(d_f-xbar(j));
%   end
%   for j=back
%     dist=x(j)-x(i)+xmax;
%     weight=exp(-dist/(n_cars*d_f/5));
%     bando_weights=bando_weights+weight;
%     xdot_bando=xdot_bando-10*weight*(d_f-xbar(j));
%   end
%   xdot(i)=xdot(i)+xdot_bando/bando_weights;
  
  xdot(i)=-10*(d_f-xbar(i));
  %xdot(i)=0.9*xdot(i);
end

plot_dynamics(x,xdot,xbar)

function plot_dynamics(x,xdot,xbar)
global tidx dt
if ~mod(tidx,500)
  figure(46)
  n_cars=numel(x)/2;
  subplot(411)
  stem(x(1:n_cars)-min(x(1:n_cars)))
  subplot(412)
  stem(xdot(1:n_cars))
  subplot(413)
  stem(xbar)
  subplot(414)
  plot_cars(x,46);
  tidx*dt
  tidx
  pause;
end

function xdot=dynamics3(x,t,u)
global xmax active_cars k_pd d_f
n_cars=numel(x)/2;
xbar1=(x(n_cars)+xmax)-x(1);
xbar=[xbar1;x(1:n_cars-1)-x(2:n_cars)];
xdot=zeros(2*n_cars,1);
xdot(1:n_cars)=x(n_cars+1:end);
xdotbar=[xdot(n_cars)-xdot(1);xdot(1:n_cars-1)-xdot(2:n_cars)];
gamma=1;
C=1;
D=1.0;
xdot(n_cars+1:end)=C*xdotbar./(0.6+xbar.^(gamma+1))-D*(d_f-xbar);

% Active cars
for i=active_cars
  if i==1
    fprintf('For now, the first car cannot be an active car\n');
  end
  xdot(i)=-k_pd*(d_f-(x(i-1)-x(i)));
end
%for i=active_cars
%  xdot(i+n_cars)=0;
%end
%xdot=min(xdot,3);

plot_dynamics(x,xdot,xbar);

function xdot=dynamics2(x,t,u)
global k_pd d_f xmax
n_cars=numel(x)/2;
xdot=zeros(2*n_cars,1);
xbar1=x(n_cars)+xmax-x(1);
xbar=[xbar1;x(1:n_cars-1)-x(2:n_cars)];
% Derivatives of position
xdot(1:n_cars)=-k_pd*(d_f-xbar);
T=25;
c=0.5;
p_dumb=0.1;
dumb_cars=logical([rand(n_cars,1)<p_dumb;zeros(n_cars,1)]);
% xdot(dumb_cars)=xdot(dumb_cars)*double(mod(t,T)>c*T);
% Double derivatives of position
xdot(n_cars+1)=-k_pd*(xdot(1)-xdot(n_cars));
xdot(n_cars+2:2*n_cars)=-k_pd*(xdot(2:n_cars)-xdot(1:n_cars-1));

plot_dynamics(x,xdot,xbar);


function x_wrap=wrap_around(x)
global xmax;
n_cars=numel(x)/2;
x_wrap=x;
for i=1:n_cars
  if x(i)>xmax
    x_wrap(i)=mod(x(i),xmax);
  end
  if x(i)<0
    x_wrap(i)=xmax-mod(-x(i),xmax);
  end
end

function plot_cars(x,fig_handle)
global xmax active_cars;
x=wrap_around(x);
n_cars=numel(x)/2;
figure(fig_handle);
scatter(x(1:n_cars),zeros(n_cars,1),'ks','SizeData',20,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
for i=active_cars
  scatter(x(i),0,'rs','SizeData',20,'MarkerEdgeColor','r','MarkerFaceColor','r');
end
hold off
xlim([0,xmax])
ylim([-0.5,0.5])
%set(gcf,'Position',[200,200,1400,800])
drawnow;

function d=wrap_diff(x1,x2)
global xmax
diffs=[x1-x2+xmax,x1-x2,x1-x2-xmax];
absdiffs=abs(diffs);
m=min(absdiffs);
mi=find(absdiffs==m);
mi=mi(1);
d=diffs(mi);
