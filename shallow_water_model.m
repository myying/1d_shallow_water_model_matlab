%clear all; close all;

%parameters
dx=100000;
dt=9900;
nx=50;
nt=20*dx/(10*dt);

pi=4*atan(1.0);
vel=10;
hbar=8000;%base level of height
amp=100;  %amplitude of hill
hc=25;    %center location of hill
hw=10;    %width of hill

g=9.8;
f=0.0001;
kdif=10;

nifcor=0;
nifwind=0;
nifdif=0;
nifad=1;

scale_factor=nx;
show_animation=1;
animation_stride=1000;
animation_delay=0.001;

%initialize
dt2=dt*2; dx2=dx*2; diffac=kdif*dt/dx/dx;
time=0;
u(1:nx)=vel;
v(1:nx)=0.0; 
h(1:nx)=hbar;
ind=ceil(hc-hw/2):floor(hc+hw/2);
h(ind)=hbar+amp*(1+cos(2*pi*(ind-hc)/hw));
%for j=jhin:jhend
%  h(j)=h(j)+amp*(1.+sin(twopi*(j-jhin)/(jhend-jhin)-0.25*twopi));
%end
h0=h;
tracer1=find(h==max(h),1); maxh(1)=max(h);
ub=u; vb=v; hb=h;
ua=u; va=v; ha=h;
x(1:nx)=(1:nx)*dx/1000;
dtcalc=dt;

%run model
figure('Position',[100 100 560 560]);
for n=1:nt
  j=2:nx-1;
  ha(j)=hb(j)-nifad*u(j).*(dtcalc/dx2).*(h(j+1)-h(j-1)) - h(j).*(dtcalc/dx2).*(u(j+1)-u(j-1));
  if(nifwind==1)
    ua(j)=ub(j)-nifad*u(j).*(dtcalc/dx2).*(u(j+1)-u(j-1)) - (g*dtcalc/dx2)*(h(j+1)-h(j-1));
  end
  if(nifdif==1)
    ha(j)=ha(j)+diffac*(hb(j+1)-2*hb(j)+hb(j-1));
    if(nifwind==1)
      ua(j)=ua(j)+diffac*(ub(j+1)-2*ub(j)+ub(j-1));
    end
  end
  if(nifcor==1)
    ua(j)=ua(j)+f*v(j)*dtcalc;
    va(j)=vb(j)-nifad*u(j).*(dtcalc/dx2).*(v(j+1)-v(j-1)) - f*u(j)*dtcalc;
  end

  %boundary condition
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
    plot(x,h); hold on;
    %analytic solution for advection equation
    x0=mod(x+vel*time/1000,max(x));
    n0=find(x0(2:end)-x0(1:end-1) <0);
    plot(x0(1:n0),h0(1:n0),'r'); 
    plot(x0(n0+1:end),h0(n0+1:end),'r'); 
    hold off; 
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
  
  if(n==floor(20*dx/vel/dt)) %trace max(h) across 20 grid points
    tracer2=find(h==max(h),1);
    Cd=(tracer2-tracer1)*dx/(floor(20*dx/vel/dt)*dt); disp(['Cd=' num2str(Cd)]);
  end
  maxh(n+1)=max(h);
  
  if(sum(any(h>1e10))>0) 
    disp(['model exploded at t=' num2str(time) 's (' num2str(n) ' steps)']); break; 
  end
end
