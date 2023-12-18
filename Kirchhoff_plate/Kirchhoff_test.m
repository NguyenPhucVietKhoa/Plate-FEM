%%
clear
close all
clc

%Dimensions
Lx=1;
Ly=1;
h=0.001;

%Material parameters
E=1.3e11;
nu=0.27;
alp=2.8e-6;

dT = 100;

N1 = 50;
N2 = 50;
lex = Lx/N1;
ley = Ly/N2;
r=lex/ley;

Ke = Stiffness_matrix_K(E,nu,lex,ley,h);
Fe = Thermal_stress_K(E,nu,alp,lex,ley,h,dT);
%number of nodes along x1:
n1=N1+1;
%number of nodes along x2:
n2=N2+1;

II = zeros(N1*N2*12^2,1);
JJ = zeros(N1*N2*12^2,1);
SSK = zeros(N1*N2*12^2,1);
index=0;

F=zeros(n1*n2*3,1);

for J=1:(n2-1)
    for I=1:(n1-1)
        DofLEe=[3*(J-1)*n1+3*(I-1)+(1:6) 3*J*n1+3*(I-1)+(4:6) 3*J*n1+3*(I-1)+(1:3)];

        F(DofLEe,1) = F(DofLEe,1) +Fe;
        for tII = 1:12
            for tJJ = 1:12
                index=index+1;
                II(index,1) = DofLEe(1,tII);
                JJ(index,1) = DofLEe(1,tJJ);
                SSK(index) = Ke(tII,tJJ);
            end
        end
    end
end

K = sparse(II,JJ,SSK,n1*n2*3,n1*n2*3);

%CLs

DOFb=zeros(n1+n2+1,1);
select_w=zeros(n1*n2,1);

DOFb(1:3)=1:3;
select_w(1)=1;
for I=2:n1
    DOFb(3+(I-1))=3*(I-1)+2;
    select_w(I)=2*(I-1);
end
for I=1:(n2-1)
    DOFb(n1+2+I)=3*n1*I+3;
    select_w(n1*I+(1:n1))=2*(n1-1)+1+(I-1)*(3*(n1-1)+2)+[1,2+(1:3:(3*n1-3))];
end

DOFg=1:(3*n1*n2);
DOFI=setdiff(DOFg,DOFb);

Ku=K(DOFI,DOFI);
Fu=F(DOFI,1);

u=Ku\Fu;
u_g=[0;u];

x_plot=0:lex:Lx;
y_plot=0:ley:Ly;

[XX,YY]=meshgrid(x_plot,y_plot);
ZZ=zeros(size(XX));

figure(1);

q_plot=zeros(n1,n2);
q_plot(1:n1,1:n2)=reshape(-u_g(select_w),n1,n2);

surf(XX,YY,ZZ)
hold on
% surf(-XX,YY,ZZ)
% surf(-XX,-YY,ZZ)
% surf(XX,-YY,ZZ)
surf(XX,YY,q_plot)
% surf(-X,Y,q_plot)
% surf(X,-Y,q_plot)
% surf(-X,-Y,q_plot)
axis([0 Lx 0 Ly])
title("Deplacement |Uz|(m) - Kirchhoff plate")
colorbar
colormap turbo
shading interp