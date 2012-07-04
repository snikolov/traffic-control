%test_traffic
function test_traffic_smoothing
global k_pd d_f xmax active_cars dt tidx tau n_cars t_h v0 n_steps skip_steps k1 k2 ka1 ka2 PLOT K dd vd A B

close all
rng('default');

PLOT=1;

k_pd=1;

% Forward simulate dynamics.
dt=0.001;

% Headway time
t_h=1;

% Control gains
k1=1;
k2=1;

ka1=0.00125;
ka2=0.00125;

n_steps=20000;
skip_steps=5;

n_cars=22;
active_cars=[];%2:2:n_cars;

if any(active_cars==1)
  error('1 cannot be an active car for now');
end
    
T=0:dt:dt*n_steps;

v0=1;
d_th=v0*t_h;
d_f=1e-1;
d_init=1*d_th;

% Positions.
% x=[linspace(n_cars*d_init,d_init,n_cars)';[1;zeros(n_cars-1,1)]];
% x=2*n_cars*d_init*rand(n_cars,1);
% x=sort(x,'descend');
x=linspace(1+(n_cars-1)*d_init,1,n_cars)';

% Velocities.
x=[x;v0;zeros(n_cars-1,1)];

% Accelerations
x=[x;zeros(n_cars,1)];

%x=[x;0.4+0.05*rand(n_cars,1)];
xmax=max(x(1:n_cars));

% Reaction Delay
delay=1;
tau=0.4;

C1=diag(-1*ones(n_cars,1))+diag(ones(n_cars-1,1),-1);
C1(1,n_cars)=1;
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
% A(1,:)=[];
% A(:,1)=[];
% A(n_cars+1,:)=[];
% A(:,n_cars+1)=[];
% A(:,2*n_cars+1)=[];
% A(2*n_cars+1,:)=[];
% active_cars=active_cars-1;
n_active=numel(active_cars);
if n_active>0
  if ~delay
    A(active_cars+n_cars,:)=0;
    B=eye(2*n_cars);
    B=B(:,n_cars+active_cars);
    Q=eye(2*n_cars);
    Q(1:n_cars+1,1:n_cars+1)=eye(n_cars+1);
  else
    A(active_cars+2*n_cars,1:2*n_cars)=0;
    B=eye(3*n_cars)/tau;
    B=B(:,2*n_cars+active_cars);
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

% Desired distances and velocities for lqr.
dd=1;
vd=1;

% sgd(x);
score=run(x)

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function sgd(x)
global k1 k2 ka1 ka2
% Stochastic Gradient Descent
n_trials=200;
PLOT=1;
scores=zeros(n_trials,1);
for trial=1:n_trials
  % Evaluate objective at current point
  score_old=run(x);
  ka1_old=ka1;
  ka2_old=ka2;
  fprintf('Current point: ka1=%.5f,ka2=%.5f,score=%.4f\n',ka1,ka2,score_old);
  % Generate test point.
  sigma=0.002; 
  for attempts=1:10
    dk=normrnd(0,sigma,2,1);
    ka1=max(ka1_old+dk(1),0);
    ka2=max(ka2_old+dk(2),0);
    % Evaluate objective at test point.
    score=run(x);
    if ~isnan(score)
      break;
    end
  end
  scores(trial)=score;
  if isnan(score)
    fprintf('No valid gain found near ka1=%.5f,ka2=%.5f\n',ka1_old,ka2_old);
  end
  fprintf('Test point: ka1=%.5f,ka2=%.5f,score=%.4f\n',ka1,ka2,score);
  figure(346)
  subplot(211)
  hold on
  scatter(ka1,ka2,'ks','SizeData',5)
  xlim([0,0.1])
  ylim([0,0.1])
  subplot(212)
  plot(scores);
  
  % Update point.
  eta=0.005;
  ka1=max(ka1_old-eta*(ka1-ka1_old)*(score-score_old),0);
  ka2=max(ka2_old-eta*(ka2-ka2_old)*(score-score_old),0);
  
  subplot(211)
  scatter(ka1,ka2)
end

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function stable=is_stable_gain(k1,k2)
global t_h tau
stable=(k2+t_h*k1<=1/(2*tau)&2*t_h*k2+t_h^2*k1>2)| ...
       (k2+t_h*k1>=1/(2*tau)&((k2-1/(2*tau))^2<(t_h/tau-2)*k1)); 
  
%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function penalty=run(x)
global dt n_steps xmax d_f n_cars PLOT
penalty=0;
T=0:dt:dt*n_steps;
global costs
cost_skip_steps=100;
% costs=zeros(ceil(numel(T)/cost_skip_steps),1);
NO_COLLISIONS=0;
NO_BACKWARD=1;
for tidx=1:numel(T)
  t=T(tidx);
  
  if PLOT
    plot_dynamics(x(1:n_cars),x(n_cars+1:2*n_cars),x(2*n_cars+1:3*n_cars),tidx);
  end
  
%   u=control_pd(x,t);
%   xdot=dynamics(x,t,u);
  d=x_to_d(x);
  u=control_lqr(d);
  ddot=dynamics_lqr(d,u);
  xdot=ddot_to_xdot(ddot);
  
  if NO_BACKWARD
    % Threshold qdot to be positive
    xdot(1:n_cars)=max(0,xdot(1:n_cars));
  end
  
  x(n_cars+1:2*n_cars)=x(n_cars+1:2*n_cars)+dt*xdot(n_cars+1:2*n_cars); 
  x(2*n_cars+1:end)=x(2*n_cars+1:end)+dt*xdot(2*n_cars+1:end);
  
  x(1)=x(1)+dt*xdot(1);
  if NO_COLLISIONS
    x(1)=min(x(1),x(n_cars)+xmax-d_f);
    if x(1)==x(n_cars)+xmax-d_f
      x(1+n_cars)=0;
      x(1+2*n_cars)=0;
    end
  end
  for i=2:n_cars
    x(i)=x(i)+dt*xdot(i);
    if NO_COLLISIONS
      x(i)=min(x(i),x(i-1)-d_f);
      % If car got stopped, set velocity and acceleration to 0.
      if x(i)==x(i-1)-d_f
        x(i+n_cars)=0;
        x(i+2*n_cars)=0;
      end
    end    
    %x(i)=x(i)+dt*xdot(i);
  end
  c=cost(x,u);
  % Only take latter half
  if tidx>2*numel(T)/3
    penalty=penalty+c;
  end
  %if ~mod(tidx,cost_skip_steps)
  %  costs(ceil(tidx/cost_skip_steps))=c;
  %end
end

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function d=x_to_d(x)
global n_cars dd vd xmax
d=zeros(3*n_cars-3,1);
d(1:n_cars-1)=x(2:n_cars)-x(1:n_cars-1)-dd;
d(n_cars:2*n_cars-2)=x(n_cars+2:2*n_cars)-x(n_cars+1:2*n_cars-1)-vd;
d(2*n_cars-1:3*n_cars-3)=x(2*n_cars+2:3*n_cars);

dtmp=zeros(n_cars*3,1);

dtmp(2:n_cars)=d(1:n_cars-1);
dtmp(n_cars+2:2*n_cars)=d(n_cars:2*n_cars-2);
dtmp(2*n_cars+2:3*n_cars)=d(2*n_cars-1:3*n_cars-3);

dtmp(1)=x(n_cars)-x(1)+xmax-dd;
dtmp(n_cars+1)=x(2*n_cars)-x(n_cars+1)-vd;
dtmp(2*n_cars-2)=x(2*n_cars+1);
d=dtmp;

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function xdot=ddot_to_xdot(ddot)
global n_cars dd vd
xdot=zeros(3*n_cars,1);
xdot(1)=0;
%xdot(2:n_cars)=cumsum(ddot(1:n_cars-1)+dd);
xdot(2:n_cars)=cumsum(ddot(2:n_cars)+dd);
xdot(n_cars+1)=0;
%xdot(n_cars+2:2*n_cars)=cumsum(ddot(n_cars:2*n_cars-2)+vd);
xdot(n_cars+2:2*n_cars)=cumsum(ddot(n_cars:2*n_cars-2)+vd);
xdot(2*n_cars+1)=0;
xdot(2*n_cars+2:3*n_cars)=cumsum(ddot(2*n_cars-1:3*n_cars-3));

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function g=cost(x,u)
global xmax t_h n_cars d_f
q=x(1:n_cars);
qbar=[q(end)+xmax-q(1);q(1:n_cars-1)-q(2:end)];
qdot=x(n_cars+1:2*n_cars);
qdotbar=[qdot(end)-q(1);qdot(1:n_cars-1)-qdot(2:end)];
R=1;
Q1=1;
Q2=0.0001;
lambda=0.01;
% Remove highest error. This error is the gap in the string of cars
% g=(Q1*norm(qbar-t_h*qdot)^2 + Q2*norm(qdotbar)^2 + R*norm(u)^2)/mean(qdot);
qbar_max=max(qbar);
qbar_max_idx=find(qbar==qbar_max);
qbar_max_idx=qbar_max_idx(1);
qbar0=qbar;
qdot0=qdot;
qbar0(qbar_max_idx)=[];
qdot0(qbar_max_idx)=[];
g=(Q1*norm(qbar0-t_h*qdot0-d_f)^2 + Q2*norm(qdotbar)^2 + R*norm(u))/(mean(qdot)+lambda);

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function u=control_lqr(d)
global K
u=-K*d;

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function ddot=dynamics_lqr(d,u)
global A B
ddot=A*d+B*u;

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function u=control_pd(x,t)
global xmax t_h n_cars d_f active_cars k1 k2 ka1 ka2
q=x(1:n_cars);
qbar1=(q(end)+xmax)-q(1);
qbar=[qbar1;q(1:end-1)-q(2:end)];
qdot=x(n_cars+1:2*n_cars);
qdotbar1=qdot(end)-qdot(1);
qdotbar=[qdotbar1;qdot(1:end-1)-qdot(2:end)];
u=k1*(qbar-t_h*qdot-d_f)+k2*qdotbar;
%u=k1*(qbar-d_f)+k2*qdotbar;
u(active_cars)=ka1*(qbar(active_cars)-t_h*qdot(active_cars)-d_f)+ka2*qdotbar(active_cars);
%u(active_cars)=ka1*(qbar(active_cars)-d_f)+ka2*qdotbar(active_cars);

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function xdot=dynamics(x,t,u)
global xmax active_cars k_pd d_f dt tau n_cars PLOT
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

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function plot_dynamics(q,qdot,qddot,tidx)
global dt d_f skip_steps costs xmax
if ~mod(tidx-1,skip_steps)
  % figure(46)
  
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
  stem([q(end)-q(1)+xmax;q(1:end-1)-q(2:end)])
  title('R')
  
  %subplot(615)
  %plot(costs)
  %title('Instantaneous cost');
  
  subplot(515)
  plot_cars(q,46);
  title(sprintf('Time: %.3f',tidx*dt));
  %title('k_1=%.3f, k_2=%.3f, k_{a,1}=%.5f, k_{a,2}=%.5f, Time: %.3f',k1,k2,ka1,ka2,tidx*dt);

  fprintf('Time: %.4f\n',tidx*dt);
  fprintf('Index: %d\n',tidx);
  
  %drawnow;
  %pause(0.01);
  %pause;
  
  %set(gcf,'Position',[200 200 1000 700])
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
global xmax active_cars n_cars;
q=wrap_around(q);
figure(fig_handle);
scatter(q,zeros(size(q)),linspace(10,80,n_cars),linspace(1,32,n_cars),'filled');
%scatter(q,zeros(size(q)),'ks','SizeData',20,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
for i=active_cars
  scatter(q(i),0,'ro','SizeData',20,'MarkerEdgeColor','r','MarkerFaceColor','r');
end
hold off
xlim([0,xmax])
ylim([-0.5,0.5])
%set(gcf,'Position',[200,200,1400,800])
colormap(cool);
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
