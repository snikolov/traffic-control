%test_traffic
function test_traffic_smoothing
global k_pd d_f xmax active_cars dt tidx tau n_cars t_h

close all

k_pd=1;
d_f=1e-1;

% Forward simulate dynamics.
dt=0.005;
% Servo loop time lag
tau=0.55;
% Headway time
t_h=1;

% Control gains
k1=2;
k2=2;

ka1=0.00125;
ka2=0.00125;

n_steps=500000;
global skip_steps
skip_steps=5000;
n_cars=50;
active_cars=[2,6,10,14,16,20,24];
T=0:dt:dt*n_steps;
v0=1;
d_th=v0*t_h;
d_init=1*d_th;

% Positions.
% x=[linspace(n_cars*d_init,d_init,n_cars)';[1;zeros(n_cars-1,1)]];
%x=.2*n_cars*d_init*rand(n_cars,1);
%x=sort(x,'descend');
x=linspace(0.001+(n_cars-1)*d_init,0.001,n_cars)';

% Velocities.
x=[x;v0;zeros(n_cars-1,1)];

% Accelerations
x=[x;zeros(n_cars,1)];

%x=[x;0.4+0.05*rand(n_cars,1)];
xmax=max(x(1:n_cars));

for tidx=1:numel(T)
  t=T(tidx);
  u=control_pd(x,t,k1,k2,ka1,ka2);
  xdot=dynamics(x,t,u);
  
  % Threshold qdot to be positive
  xdot(1:n_cars)=max(0,xdot(1:n_cars));
  
  x(1)=min(x(1)+dt*xdot(1),x(n_cars)+xmax-d_f);
  for i=2:n_cars
    x(i)=min(x(i)+dt*xdot(i),x(i-1)-d_f);
    %x(i)=x(i)+dt*xdot(i);
  end
  x(n_cars+1:2*n_cars)=x(n_cars+1:2*n_cars)+dt*xdot(n_cars+1:2*n_cars);
  x(2*n_cars+1:end)=x(2*n_cars+1:end)+dt*xdot(2*n_cars+1:end);
end

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function u=control_pd(x,t,k1,k2,ka1,ka2)
global xmax t_h n_cars d_f active_cars
q=x(1:n_cars);
qbar1=(q(end)+xmax)-q(1);
qbar=[qbar1;q(1:end-1)-q(2:end)];
qdot=x(n_cars+1:2*n_cars);
qdotbar1=qdot(end)-qdot(1);
qdotbar=[qdotbar1;qdot(1:end-1)-qdot(2:end)];
u=k1*(qbar-t_h*qdot)+k2*qdotbar;
%u=k1*(qbar-d_f)+k2*qdotbar;
u(active_cars)=ka1*(qbar(active_cars)-t_h*qdot(active_cars))+ka2*qdotbar(active_cars);
%u(active_cars)=ka1*(qbar(active_cars)-d_f)+ka2*qdotbar(active_cars);

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function xdot=dynamics(x,t,u)
global xmax active_cars k_pd d_f dt xdothist tau n_cars
q=x(1:n_cars);
qdot=x(n_cars+1:2*n_cars);
qddot=x(2*n_cars+1:3*n_cars);
xdot=[qdot;qddot;(u-qddot)/tau];

% % Dumb cars
% T=50*t;
% c=0.5;
% dumb_cars=logical(zeros(n_cars,1));
% dumb_cars(10)=1;
% %modulator=double(mod(t,T)>c*T);
% modulator=(1+sin(2*pi*t/T))/2;
% xdot(dumb_cars)=xdot(dumb_cars)*modulator;

plot_dynamics(q,qdot,qddot)

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function plot_dynamics(q,qdot,qddot)
global tidx dt d_f skip_steps
if ~mod(tidx-1,skip_steps)
  figure(46)
  
  subplot(511)
  stem(q-min(q))
  title('x')
  
  subplot(512)
  stem(qdot)
  title('v')
  
  subplot(513)
  stem(qddot)
  title('a')
  
  subplot(514)
  stem([q(end)-q(1);q(1:end-1)-q(2:end)])
  title('R')
  
  subplot(515)
  plot_cars(q,46);
  title(sprintf('Time: %.3f',tidx*dt));
  %title('k_1=%.3f, k_2=%.3f, k_{a,1}=%.5f, k_{a,2}=%.5f, Time: %.3f',k1,k2,ka1,ka2,tidx*dt);

  fprintf('Time: %.4f\n',tidx*dt);
  fprintf('Index: %d\n',tidx);
  pause(0.01);
  pause;
end

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function q_wrap=wrap_around(q)
global xmax;
q_wrap=q;
for i=1:numel(q)
  if q(i)>xmax
    q_wrap(i)=mod(q(i),xmax);
  end
  if q(i)<0
    q_wrap(i)=xmax-mod(-q(i),xmax);
  end
end

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function plot_cars(q,fig_handle)
global xmax active_cars;
q=wrap_around(q);
figure(fig_handle);
scatter(q,zeros(size(q)),'ks','SizeData',20,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
for i=active_cars
  scatter(q(i),0,'rs','SizeData',20,'MarkerEdgeColor','r','MarkerFaceColor','r');
end
hold off
xlim([0,xmax])
ylim([-0.5,0.5])
%set(gcf,'Position',[200,200,1400,800])
drawnow;

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function d=wrap_diff(x1,x2)
global xmax
diffs=[x1-x2+xmax,x1-x2,x1-x2-xmax];
absdiffs=abs(diffs);
m=min(absdiffs);
mi=find(absdiffs==m);
mi=mi(1);
d=diffs(mi);

% LocalWords:  linspace init rand
