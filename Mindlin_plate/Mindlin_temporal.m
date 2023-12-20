clear
%close all
clc

%Dimensions
L1=1;
L2=1;
L3=0.001;

%Material parameters
E=2.1e11;
rho=7800;
nu=0.3;
k = 5/6;
%numbers of elements
N1=60;
N2=60;

%Element stiffness matrix
l1=L1/N1;
l2=L2/N2;
r=l1/l2;

Ke = Stiffness_matrix_M(E,nu,k,l1,l2,L3);

%Element mass matrix

Me = Mass_matrix_M(rho,l1,l2,L3);
% Assembly
ntot = (2*N1+1)*(N2+1)+(N1+1)*N2;

II = zeros(N1*N2*24^2,1);
JJ = zeros(N1*N2*24^2,1);
SSK = zeros(N1*N2*24^2,1);
SSM = zeros(N1*N2*24^2,1);
index=0;
for J = 1:N2
    for I = 1:N1
        DofLEe = 6*(I-1)+3*(3*N1+2)*(J-1)+ ...
                [(1:3) (7:9) 3*(3*N1+2)+(7:9) 3*(3*N1+2)+(1:3) ...
                 (4:6) 3*(2*N1+1)-3*(I-1)+(4:6) 3*(3*N1+2)+(4:6) 3*(2*N1+1)-3*(I-1)+(1:3)];
        for tII = 1:24
            for tJJ = 1:24
                index=index+1;
                II(index,1) = DofLEe(1,tII);
                JJ(index,1) = DofLEe(1,tJJ);
                SSK(index) = Ke(tII,tJJ);
                SSM(index) = Me(tII,tJJ);
            end
        end
    end
end

Kgs = sparse(II,JJ,SSK,ntot*3,ntot*3);
Mgs = sparse(II,JJ,SSM,ntot*3,ntot*3);

%kinematic constraints
n1 = 2*N1+1;
n2 = 2*N2+1;
DOFb=zeros(2*n1+2*(n2-2),1);
select_w=zeros(n1*n2-N1*N2-(2*n1+2*(n2-2)),1);

for I=1:n1
    DOFb(I)=3*(I-1)+1;
    DOFb(n1+2*(n2-2)+I)=3*(n1*n2-N1*N2-n1)+3*(I-1)+1;
end
for J = 1:N2-1
    DOFb(n1+4*(J-1)+(1:2)) = 3*n1*J+3*(N1+1)*(J-1) +[1 3*(N1+1)+1];
    DOFb(n1+4*(J-1)+(3:4)) = 3*n1*J+3*(N1+1)*J-3 +[1 3*n1+1];
    select_w((N1+n1-3)*(J-1)+(1:(N1+n1-3))) = ...
       2*n1+2+(3*(N1+n1+1)-4)*(J-1)+[(3:3:(3*(N1+1)-6)) 3*(N1+1)+(1:3:(3*(n1-2)-2))];
end
DOFg=1:(3*n1*n2-3*N1*N2);
DOFI=setdiff(DOFg,DOFb);

M=Mgs(DOFI,DOFI);
K=Kgs(DOFI,DOFI);
C=0.001*M+2e-7*K;

xex=30*l1;
yex=30*l2;

Iex=1+xex/l1;
Jex=1+yex/l2;
Dofex=round(2*n1+3*(N1-1)+4+(Jex-1)*(8+3*(n1-2+N1-1))+5+(Iex-1)*6+1);


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

x1_plot=0:l1/2:L1;
x2_plot=0:l2/2:L2;

[~,N,~] = Shape_function_M(0,0);
Nr=[N(1) N(5) N(2) N(8) N(6) N(4) N(7) N(3)];

utot = zeros(3*ntot,1);
w_select = zeros(ntot,1);
u_x=zeros(N1*N2,1);
w_plot=zeros(n1*n2,1);

q_plot=zeros(n1,n2);
h = surf(x2_plot,x1_plot,real(q_plot),'facecolor','interp','edgecolor','none');
axis([0 L2 0 L1 -1e-5 1e-5]);
set(gca,'color','none')
colormap turbo
axis off
gif('Mindlin.gif')

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

    utot(DOFI,1) = u;
    %figure(1)
    for I = 1:ntot
        w_select(I,1) = utot(3*(I-1)+1,1);
    end
    a=0;
    for J=1:N2
        for I=1:N1
            a=a+1;
            Dof = (3*N1+2)*(J-1)+ ...
                  [2*(I-1)+(1:3) (I-1)+(2*N1+1)+(1:2) 2*(I-1)+(3*N1+2)+(1:3)];
            u_x(a)=Nr*w_select(Dof); 
        end
    end
    a=0;
    Do=zeros(N1*N2,1);
    for J=1:N2
        for I=1:N1
            a=a+1;
            Do(a) = 2*(I-1)+(2*N1+1)*2*(J-1)+(2*N1+3);
        end
    end
    w_plot(Do,1)=u_x;
    Di=setdiff(1:n1*n2,Do);
    w_plot(Di,1)=w_select;

    q_plot=zeros(n1,n2);
    q_plot(1:n1,1:n2)=reshape(w_plot,n1,n2);
    set(h,'Zdata',real(q_plot));
    gif
end
