%test_traffic
function test_traffic_smoothing
global k_pd d_f xmax active_cars dt tidx tau n_cars t_h v0 n_steps skip_steps k1 k2 ka1 ka2 PLOT

close all

PLOT=1;

k_pd=1;

% Forward simulate dynamics.
dt=0.1;
% Servo loop time lag
tau=0.55;
% Headway time
t_h=1;

% Control gains
k1=1;
k2=1;

ka1=0.00;
ka2=0.00;

n_steps=1000;
skip_steps=100;
n_cars=50;
active_cars=[3];%[2,6,10,14,16,20,24];
T=0:dt:dt*n_steps;
v0=1;
d_th=v0*t_h;
d_f=1;
d_init=1*d_th;

% Positions.
% x=[linspace(n_cars*d_init,d_init,n_cars)';[1;zeros(n_cars-1,1)]];
x=2*n_cars*d_init*rand(n_cars,1);
x=sort(x,'descend');
% x=linspace(10+(n_cars-1)*d_init,10,n_cars)';

% Velocities.
x=[x;v0;zeros(n_cars-1,1)];

% Accelerations
x=[x;zeros(n_cars,1)];

%x=[x;0.4+0.05*rand(n_cars,1)];
xmax=max(x(1:n_cars));

% sgd(x);
score=run(x)

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function sgd(x)
global k1 k2 ka1 ka2
ka1=k1;
ka2=k2;
% Stochastic Gradient Descent
n_trials=200;
PLOT=1;
scores=zeros(n_trials,1);
for trial=1:n_trials
  % Evaluate objective at current point
  score_old=run(x);
  ka1_old=ka1;
  ka2_old=ka2;
  % Generate test point.
  sigma=0.1; 
  for attempts=1:10
    ka1=ka1_old+normrnd(0,sigma);
    ka2=ka2_old+normrnd(0,sigma);
    ka1=max(0,ka1);
    ka2=max(0,ka2);
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
  fprintf('ka1=%.5f,ka2=%.5f,score=%.4f\n',ka1,ka2,score);
  figure(346)
  subplot(211)
  hold on
  scatter(ka1,ka2)
  xlim([0,1.5])
  ylim([0,1.5])
  subplot(212)
  plot(scores);
  % Update point.
  eta=0.0001;
  ka1=ka1_old-eta*(ka1-ka1_old)*(score-score_old);
  ka2=ka2_old-eta*(ka2-ka2_old)*(score-score_old);
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
for tidx=1:numel(T)
  t=T(tidx);
  
  if PLOT
    plot_dynamics(x(1:n_cars),x(n_cars+1:2*n_cars),x(2*n_cars+1:3*n_cars),tidx);
  end
  
  u=control_pd(x,t);
  xdot=dynamics(x,t,u);
  
  % Threshold qdot to be positive
  xdot(1:n_cars)=max(0,xdot(1:n_cars));
  
  x(n_cars+1:2*n_cars)=x(n_cars+1:2*n_cars)+dt*xdot(n_cars+1:2*n_cars); 
  x(2*n_cars+1:end)=x(2*n_cars+1:end)+dt*xdot(2*n_cars+1:end);
  
  x(1)=min(x(1)+dt*xdot(1),x(n_cars)+xmax-d_f);
  if x(1)==x(n_cars)+xmax-d_f
    x(1+n_cars)=0;
    x(1+2*n_cars)=0;
  end
  for i=2:n_cars
    x(i)=min(x(i)+dt*xdot(i),x(i-1)-d_f);
    % If car got stopped, set velocity and acceleration to 0.
    if x(i)==x(i-1)-d_f
      x(i+n_cars)=0;
      x(i+2*n_cars)=0;
    end
    
    %x(i)=x(i)+dt*xdot(i);
  end
  penalty=penalty+cost(x,u);
end

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function g=cost(x,u)
global xmax t_h n_cars
q=x(1:n_cars);
qbar=[q(end)+xmax-q(1);q(1:n_cars-1)-q(2:end)];
qdot=x(n_cars+1:2*n_cars);
qdotbar=[qdot(end)-q(1);qdot(1:n_cars-1)-qdot(2:end)];
R=1;
Q1=1;
Q2=0.0001;
% Remove highest error. This error is the gap in the string of cars
g=(Q1*norm(qbar-t_h*qdot)^2 + Q2*norm(qdotbar)^2 + R*norm(u)^2)/mean(qdot);

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
global dt d_f skip_steps
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
  %pause;
  
  set(gcf,'Position',[200 200 1000 500])
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
