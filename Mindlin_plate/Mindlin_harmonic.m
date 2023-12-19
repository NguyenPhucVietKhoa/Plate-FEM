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
k = 5/6;
%numbers of elements
N1=40;
N2=40;

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
for I=1:(n2-2)
    DOFb(n1+(I-1)*2+1)=3*n1+3*(N1+1)*(I-1)+1;
    DOFb(n1+(I-1)*2+2)=3*n1*I+3*(N1+1)*I-3+1;
end

for J = 1:N2-1
    select_w((N1+n1-3)*(J-1)+(1:(N1+n1-3))) = ...
       2*n1+2+(3*(N1+n1+1)-4)*(J-1)+[(3:3:(3*(N1+1)-6)) 3*(N1+1)+(1:3:(3*(n1-2)-2))];
end
DOFg=1:(3*n1*n2-3*N1*N2);
DOFI=setdiff(DOFg,DOFb);

M=Mgs(DOFI,DOFI);
K=Kgs(DOFI,DOFI);
C=0.001*M;

xex=5*l1;
yex=5*l2;

Iex=1+xex/l1;
Jex=1+yex/l2;
Dofex=round(2*n1+3*(N1-1)+4+(Jex-1)*(8+3*(n1-2+N1-1))+5+(Iex-1)*6+1);

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