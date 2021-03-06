% Optimal Velocity Traffic Model

function test_traffic_self_sustaining
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
a=1.8;%0.786*(2-0.1)

% rng('default');

vmax=1;
n_cars=100;
L=2.5*n_cars;

pos=flipud(linspace(L/n_cars,L,n_cars)');

pert=0.05*L/n_cars*(2*rand(n_cars,1)-1);
%pert=0.3*L/n_cars*(rand(n_cars,1)<0.01);
pos=pos+pert;

vel=vopt(xtod(pos));
%pert=0.3*L/n_cars*(rand(n_cars,1)<0.01);
vel=vel+pert;
x=[pos;vel];

%x=[pos;zeros(n_cars,1)];
%x(floor(1.1*n_cars))=1;
%x(floor(1.5*n_cars))=1;
%x(floor(1.8*n_cars))=1;

%car_ind=2:n_cars;
%active=car_ind(rand(1,n_cars-1)<0.1);
%active=[1:n_cars];
%active=[1:5:n_cars];
active=[1:25:n_cars]
%active=[];

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function run(x)
global n_cars tidx fftmag L
% Boundary for plotting
dt=0.1;
figure
iter=500000;
skip=5000;
time_evol=zeros(n_cars*floor(iter/skip),2);
fftmag=[];
for tidx=1:iter
  xdot=dynamics(x);
  if ~mod(tidx-1,skip)
    plot_cars(x,xdot);
    
    %time_evol(n_cars*(tidx/skip-1)+1:n_cars*(tidx/skip),:)=horzcat(mod(x(1:n_cars),L),tidx*dt*ones(n_cars,1));
    fftmag=[fftmag;norm(fft(xtod(x(1:n_cars))-L/n_cars))];
    
    % c=h^0.25, 1:5:end
    %cfss=ones(n_cars,1)*1.436237265147213;
    %cfss(1:5:end)=4.255050939411025;
    % c=h^0.5, 1:5:end
    %cfss=ones(n_cars,1)*1.741657386771885;
    %cfss(1:5:end)=   3.033370452897072;
    % c=h^0.5, 1:25:end
    %cfss=ones(n_cars,1)*1.928388277163331;
    %cfss(1:25:end)=3.718681347500960;
    % c=h^0.25, 1:25:end
    %cfss=ones(n_cars,1)*1.719275374745423;
    %cfss(1:25:end)=8.737391006084181;
    %fftmag=[fftmag;norm(fft(xtod(x(1:n_cars))-cfss))];
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

spatiotemporal=0;
if spatiotemporal
  scatter(time_evol(:,1),time_evol(:,2),'ks','SizeData',1);
  title(sprintf('L=%d, N=%d, a=%.4f',L,n_cars,a));
  xlabel('Position')
  ylabel('Time')
  set(gcf,'Position',[200,200,300,240]);
end
save('vmctrl_k25_fft_25.mat','fftmag')

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

xdot(n_cars+active)=a*(vopt(d(active).^0.25)-x(n_cars+active));

%vbar=[x(2*n_cars)-x(n_cars+1);x(n_cars+1:2*n_cars-1)-x(n_cars+2:2*n_cars)];
%xdot(n_cars+active)=a*(vopt(d(active))-x(n_cars+active))+1*vbar(active);

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function d=xtod(x)
global L n_cars
d=[x(n_cars)-x(1)+L;x(1:n_cars-1)-x(2:n_cars)];

%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
function v=vopt(h)
global vmax
%v=zeros(size(h));
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
global n_cars L tidx active fftmag

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
semilogy(fftmag)

% NFFT=2^nextpow2(n_cars);
% X=fft(xtod(x(1:n_cars))-L/n_cars,NFFT)/n_cars;
% freq=2*pi/2*linspace(0,1,NFFT/2+1);
% %v=x(1+n_cars:2*n_cars);
% %absfft=abs(fft(v-mean(v)));
% plot(freq,2*abs(X(1:NFFT/2+1)),'k')
% hold on

%plot(xdot(n_cars+1:2*n_cars),'k')
%title('Accelerations')

drawnow;
set(gcf,'Position',[100,100,1000,500]);
pause(0.01)
