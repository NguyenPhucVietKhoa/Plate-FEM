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
N1=40;
N2=40;

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
C=0.001*M;

xex=5*l1;
yex=5*l2;

% xmes=0.4;
% ymes=0.4;

Iex=1+xex/l1;
Jex=1+yex/l2;
Dofex=round(2*n1+(Jex-1)*(2*2+(n1-2)*3)+2+(Iex-1)*3+1);


f=0.1:0.1:30;
w=2*pi*f;

F=zeros(length(DOFI),1);
F(Dofex)=0.3;

u_mes=zeros(length(f),1);

for o=1:length(f)
    o
    D=-w(o)^2*M+1i*w(o)*C+K;
    u=D\F;
    u_mes(o)=u(Dofex);
end

figure(1)
semilogy(f,abs(u_mes));


x1_plot=0:l1:L1;
x2_plot=0:l2:L2;
%--------------------------------------------------------------------------
figure(2);
f=41.9;
w=2*pi*f;
D=-w^2*M+1i*w*C+K;
u=D\F;

q_plot=zeros(n1,n2);
q_plot(2:(n1-1),2:(n2-1))=reshape(u(select_w),(n1-2),(n2-2));
          
surf(x2_plot,x1_plot,real(q_plot),'facecolor','interp','edgecolor','none')
axis([0 L2 0 L1 -5e-3 5e-3])
title('41.9 Hz')











