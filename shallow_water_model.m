clear all; close all;

%parameters
dx=100000;  %grid spacing (m)
dt=1000;    %time step (s)
nx=50;      %number of grid points
nt=1000;    %number of time steps

pi=4*atan(1.0);
vel=10;     %advection velocity
hbar=8000;  %base level of height
amp=100;    %amplitude of hill
hc=25;      %center location of hill
hw=10;      %width of hill

g=9.8;      %gravity
f=0.0001;   %Coriolis
kdif=10;    %diffuction coef

nifcor=0;   %if 1, include Coriolis
nifwind=0;  %if 1, include other terms besides Coriolis in u,v eqns.
nifdif=0;   %if 1, include diffusion
nifad=1;    %if 1, include advection

scale_factor=nx;       %scaling factor for plotting wind vectors
show_animation=1;      %if 1, show animation
animation_stride=1;    %if =n, skip every n frames
animation_delay=0.001; %deley in between frame in seconds

%initialize
dt2=dt*2; dx2=dx*2; diffac=kdif*dt/dx/dx;
time=0;  %time in seconds
u(1:nx)=vel;
v(1:nx)=0.0; 
h(1:nx)=hbar;
ind=ceil(hc-hw/2):floor(hc+hw/2);
h(ind)=hbar+amp*(1+cos(2*pi*(ind-hc)/hw));
ub=u; vb=v; hb=h;
ua=u; va=v; ha=h;
h0=h;
loc1=find(h==max(h),1); maxh(1)=max(h);
x(1:nx)=(1:nx)*dx/1000;
dtcalc=dt;

%run model
figure('Position',[100 100 560 560]);
for n=1:nt
  %CTCS scheme (leapfrog)
  j=2:nx-1;
  ha(j)=hb(j)-nifad*u(j).*(dtcalc/dx2).*(h(j+1)-h(j-1)) ...
             -h(j).*(dtcalc/dx2).*(u(j+1)-u(j-1)) ...
             +nifdif*diffac*(hb(j+1)-2*hb(j)+hb(j-1));
  ua(j)=ub(j)-nifwind*nifad*u(j).*(dtcalc/dx2).*(u(j+1)-u(j-1)) ...
             -nifwind*(g*dtcalc/dx2)*(h(j+1)-h(j-1)) ...
             +nifcor*f*v(j)*dtcalc ...
             +nifwind*nifdif*diffac*(ub(j+1)-2*ub(j)+ub(j-1));
  va(j)=vb(j)-nifwind*nifad*u(j).*(dtcalc/dx2).*(v(j+1)-v(j-1)) ...
             -nifcor*f*u(j)*dtcalc ...
             +nifwind*nifdif*diffac*(vb(j+1)-2*vb(j)+vb(j-1));

  %boundary condition (cyclic)
  ua(1)=ua(nx-1); ua(nx)=ua(2);
  va(1)=va(nx-1); va(nx)=va(2);
  ha(1)=ha(nx-1); ha(nx)=ha(2);

  %update
  ub=u; u=ua;
  vb=v; v=va;
  hb=h; h=ha;
  time=time+dt;
  dtcalc=dt2;
  diffac=kdif*dtcalc/dx/dx;

  %plot results
  if(show_animation==1 && mod(time,animation_stride*dt)==0)
    subplot(2,1,1)
    plot(x,h);
    %analytic solution for advection equation
    %  hold on;
    %  x0=mod(x+vel*time/1000,max(x));
    %  n0=find(x0(2:end)-x0(1:end-1) <0);
    %  plot(x0(1:n0),h0(1:n0),'r'); 
    %  plot(x0(n0+1:end),h0(n0+1:end),'r'); hold off; 
    axis([min(x) max(x) hbar-amp hbar+2*amp]);
    ylabel('height'); title(['t = ' num2str(time) ' s'])
    xlabel('x (km)')
    subplot(2,1,2)
    i=1:ceil(nx/30):nx;
    quiver(x(i),0*x(i),u(i)*scale_factor,v(i)*scale_factor,'autoscale','off')
    axis equal; axis([min(x) max(x) -max(x)/5 max(x)/5]);
    xlabel('x (km)'); ylabel('y (km)');
    pause(animation_delay);
  end
  
  m=20;
  if(n==floor(m*dx/vel/dt)) %trace max(h) across m grid points
    loc2=find(h==max(h),1);
    Cd=(loc2-loc1)*dx/(floor(m*dx/vel/dt)*dt); disp(['Cd=' num2str(Cd)]);
  end
  maxh(n+1)=max(h); %record h peak amplitude time series
  
  if(sum(any(h>1e10))>0) 
    disp(['model exploded at t=' num2str(time) 's (' num2str(n) ' steps)']); break; 
  end
end
