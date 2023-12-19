clear
close all
clc

%Dimensions
L1=1;
L2=1;
L3=0.001;

%Material parameters
E=2.1e11;
rho=7800;
nu=0.3;

%numbers of elements
N1=60;
N2=60;

%Element stiffness matrix
l1=L1/N1;
l2=L2/N2;
r=l1/l2;

Ke = Stiffness_matrix_K(E,nu,l1,l2,L3);

%Element mass matrix

Me = Mass_matrix_K(rho,l1,l2,L3);

%number of nodes along x1:
n1=N1+1;
%number of nodes along x2:
n2=N2+1;


Mgs=sparse(n1*n2*3,n1*n2*3);
Kgs=sparse(n1*n2*3,n1*n2*3);

for J=1:(n2-1)
    for I=1:(n1-1)
        DofLEe=[3*(J-1)*n1+3*(I-1)+(1:6) 3*J*n1+3*(I-1)+(4:6) 3*J*n1+3*(I-1)+(1:3)];
        
        Mgs(DofLEe,DofLEe)=Mgs(DofLEe,DofLEe)+Me;
        Kgs(DofLEe,DofLEe)=Kgs(DofLEe,DofLEe)+Ke;
    end
end


%kinematic constraints
DOFb=zeros(2*n1+2*(n2-2),1);
select_w=zeros((n1-2)*(n2-2),1);

for I=1:n1
    DOFb(I)=3*(I-1)+1;
    DOFb(n1+2*(n2-2)+I)=n1*(n2-1)*3+3*(I-1)+1;
end
for I=1:(n2-2)
    DOFb(n1+(I-1)*2+1)=3*n1*I+1;
    DOFb(n1+(I-1)*2+2)=3*n1*(I+1)-3+1;
    select_w((n1-2)*(I-1)+(1:(n1-2)))=2*n1+(I-1)*(3*n1-2)+(3:3:(3*n1-6));
end
DOFg=1:(3*n1*n2);
DOFI=setdiff(DOFg,DOFb);


M=Mgs(DOFI,DOFI);
K=Kgs(DOFI,DOFI);
C=0.001*M+2e-7*K;

xex=30*l1;
yex=30*l2;

% xmes=0.4;
% ymes=0.4;

Iex=1+xex/l1;
Jex=1+yex/l2;
Dofex=round(2*n1+(Jex-1)*(2*2+(n1-2)*3)+2+(Iex-1)*3+1);

deltat=1e-4;
time=0:deltat:0.05;

f0=500;
t0=5e-3;
%l'ondelette de chapeau mexicain avec un retard de temps t0
F=zeros(length(DOFI),length(time));
F(Dofex,:)=1*(1-2*pi^2*f0^2*(time-t0).^2).*exp(-pi^2*f0^2.*(time-t0).^2);

%Newmark parameters
gamma=1/2;
beta=1/4;

%initial conditions
u0=zeros(length(DOFI),1);
u_dot0=zeros(length(DOFI),1);
u_ddot0=M\(F(:,1)-C*u_dot0-K*u0);

S=M+gamma*deltat*C+beta*deltat^2*K;

x1_plot=0:l1:L1;
x2_plot=0:l2:L2;

q_plot=zeros(n1,n2);
h = surf(x2_plot,x1_plot,real(q_plot),'facecolor','interp','edgecolor','none');
axis([0 L2 0 L1 -1e-5 1e-5]);
axis off
gif('Kirchhoff.gif')

for t=2:length(time)
    if t==2
        u_dot_star=u_dot0+(1-gamma)*deltat*u_ddot0;
        u_star=u0+deltat*u_dot0+(1/2-beta)*deltat^2*u_ddot0;
    else
        u_dot_star=u_dot+(1-gamma)*deltat*u_ddot;
        u_star=u+deltat*u_dot+(1/2-beta)*deltat^2*u_ddot;
    end
    u_ddot=S\(F(:,t)-C*u_dot_star-K*u_star);
    u_dot=u_dot_star+gamma*deltat*u_ddot;
    u=u_star+beta*deltat^2*u_ddot;
    
    %figure(1)
    q_plot=zeros(n1,n2);
    q_plot(2:(n1-1),2:(n2-1))=reshape(u(select_w),(n1-2),(n2-2));
    set(h,'Zdata',real(q_plot));
    gif
    %surf(x2_plot,x1_plot,real(q_plot),'facecolor','interp','edgecolor','none')
    %axis([0 L2 0 L1 -5e-5 5e-5]);
    %pause(0.001);
    
    %figure(2)
    % plot(time(1:t),F(Dofex,1:t),'LineWidth',2);
    % set(h2, 'Position', [100, 200, 1300, 600])
    % axis([0 time(length(time)) -5 5])
    %pause(0.001);
end