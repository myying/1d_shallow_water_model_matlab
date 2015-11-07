clear all; close all;

%parameters
dx=100000;  %grid spacing (m)
dt=1000;    %time step (s)
nx=52;      %number of grid points (including 2 ghost nodes)
nt=1000;    %number of time steps

pi=4*atan(1.0);
C=10;       %advection velocity
hbar=8000;  %base level of height
amp=100;    %amplitude of hill
hc=25;      %center location of hill
hw=10;      %width of hill

g=9.8;      %gravity
f=0.0001;   %Coriolis
kdif=10;    %diffuction coef

numerics=1; %1:CTCS, 2:FTCS, 3: FTBS numerical schemes

nifcor=0;   %if 1, include Coriolis
nifwind=0;  %if 1, include other terms besides Coriolis in u,v eqns.
nifdif=0;   %if 1, include diffusion
nifad=1;    %if 1, include advection

scale_factor=nx;       %scaling factor for plotting wind vectors
show_animation=1;      %if 1, show animation; turn this off when profiling.
animation_stride=10;   %if =n, skip every n frames
animation_delay=0.001; %delay in between frame in seconds
show_analytic=1;       %show analytical solution (linear advection eqn)
show_wind=1;           %show wind vector

%initial condition
u(1:nx)=C;
v(1:nx)=0.0; 
h(1:nx)=hbar;
ind=ceil(hc-hw/2):floor(hc+hw/2);
h(ind)=hbar+amp*(1+cos(2*pi*(ind-hc)/hw));

h0=h;
loc1=find(h==max(h),1); 
maxh(1:nt+1)=NaN; maxh(1)=max(h);
x(1:nx)=(1:nx)*dx/1000;

figure('Position',[100 100 560 560]);

%run model
time=0;  %time in seconds
t_start=cputime;

for n=1:nt
  j=2:nx-1;
  if(numerics==1) %CTCS scheme (leapfrog)
    if(n==1)
      ub=u; vb=v; hb=h;
      ua=u; va=v; ha=h;
      dtcalc=dt;
    else
      dtcalc=dt*2;
    end
    diffac=kdif*dtcalc/dx/dx;
    ha(j)=hb(j)-nifad*u(j).*(dtcalc/(dx*2)).*(h(j+1)-h(j-1)) ...
               -h(j).*(dtcalc/(dx*2)).*(u(j+1)-u(j-1)) ...
               +nifdif*diffac*(h(j+1)-2*h(j)+h(j-1));
    ua(j)=ub(j)-nifwind*nifad*u(j).*(dtcalc/(dx*2)).*(u(j+1)-u(j-1)) ...
               -nifwind*(g*dtcalc/(dx*2))*(h(j+1)-h(j-1)) ...
               +nifcor*f*v(j)*dtcalc ...
               +nifwind*nifdif*diffac*(u(j+1)-2*u(j)+u(j-1));
    va(j)=vb(j)-nifwind*nifad*u(j).*(dtcalc/(dx*2)).*(v(j+1)-v(j-1)) ...
               -nifcor*f*u(j)*dtcalc ...
               +nifwind*nifdif*diffac*(v(j+1)-2*v(j)+v(j-1));
    ub=u; u=ua;
    vb=v; v=va;
    hb=h; h=ha;

  elseif(numerics==2) %FTCS scheme
    diffac=kdif*dt/dx/dx;
    ha(j)=h(j)-nifad*u(j).*(dt/(dx*2)).*(h(j+1)-h(j-1)) ...
              -h(j).*(dt/(dx*2)).*(u(j+1)-u(j-1)) ...
              +nifdif*diffac*(h(j+1)-2*h(j)+h(j-1));
    ua(j)=u(j)-nifwind*nifad*u(j).*(dt/(dx*2)).*(u(j+1)-u(j-1)) ...
              -nifwind*(g*dt/(dx*2))*(h(j+1)-h(j-1)) ...
              +nifcor*f*v(j)*dt ...
              +nifwind*nifdif*diffac*(u(j+1)-2*u(j)+u(j-1));
    va(j)=v(j)-nifwind*nifad*u(j).*(dt/(dx*2)).*(v(j+1)-v(j-1)) ...
              -nifcor*f*u(j)*dt ...
              +nifwind*nifdif*diffac*(v(j+1)-2*v(j)+v(j-1));
    h=ha; u=ua; v=va;
  
  elseif(numerics==3) %FTBS scheme
    diffac=kdif*dt/dx/dx;
    ha(j)=h(j)-nifad*u(j).*(dt/dx).*(h(j)-h(j-1)) ...
              -h(j).*(dt/dx).*(u(j)-u(j-1)) ...
              +nifdif*diffac*(h(j+1)-2*h(j)+h(j-1));
    ua(j)=u(j)-nifwind*nifad*u(j).*(dt/dx).*(u(j)-u(j-1)) ...
              -nifwind*(g*dt/dx)*(h(j)-h(j-1)) ...
              +nifcor*f*v(j)*dt ...
              +nifwind*nifdif*diffac*(u(j+1)-2*u(j)+u(j-1));
    va(j)=v(j)-nifwind*nifad*u(j).*(dt/dx).*(v(j)-v(j-1)) ...
              -nifcor*f*u(j)*dt ...
              +nifwind*nifdif*diffac*(v(j+1)-2*v(j)+v(j-1));
    h=ha; u=ua; v=va;
  end
  
  time=time+dt;
  
  %boundary condition (cyclic)
  u(1)=u(nx-1); u(nx)=u(2);
  v(1)=v(nx-1); v(nx)=v(2);
  h(1)=h(nx-1); h(nx)=h(2);
  
  %plot results
  if(show_animation==1 && mod(n,animation_stride)==0)
    subplot(2,1,1)
    plot(x,h);
    if(show_analytic==1) %analytic solution for advection equation
      hold on;
      x0=mod(x+C*time/1000,(nx-2)*dx/1000);
      n0=find(x0(2:end-1)-x0(1:end-2)<0);
      plot(x0(1:n0),h0(1:n0),'r');
      plot(x0(n0+1:end-1),h0(n0+1:end-1),'r');
      hold off;
    end
    axis([dx/1000 (nx-2)*dx/1000 hbar-amp hbar+2*amp]);
    ylabel('height'); title(['t = ' num2str(time) ' s'])
    xlabel('x (km)')
    if(show_wind==1)
      subplot(2,1,2)
      i=1:ceil(nx/30):nx;
      quiver(x(i),0*x(i),u(i)*scale_factor,v(i)*scale_factor,'autoscale','off')
      axis equal; axis([dx (nx-2)*dx -(nx-2)*dx/5 (nx-2)*dx/5]/1000);
      xlabel('x (km)'); ylabel('y (km)');
    end
    pause(animation_delay);
  end
  
  m=20;
  if(n==floor(m*dx/C/dt)) %trace max(h) across m grid points
    loc2=find(h==max(h),1);
    Cd=(loc2-loc1)*dx/(floor(m*dx/C/dt)*dt); disp(['Cd=' num2str(Cd)]);
  end
  maxh(n+1)=max(h); %record h peak amplitude time series
  
  if(sum(any(h>1e10))>0) 
    disp(['model exploded at t=' num2str(time) 's (' num2str(n) ' steps)']); break; 
  end
end

t_end=cputime;
disp(['model integration takes ' num2str(t_end-t_start) ' s']);
