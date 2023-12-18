%%
clear
close all
clc

%Dimensions
Lx=0.01;
Ly=0.01;
h=500e-6;

%Material parameters
E=1.3e11;
nu=0.27;
alp=2.8e-6;
k = 5/6;

dT = 100;

N1 = 50;
N2 = 50;
lex = Lx/N1;
ley = Ly/N2;

%%
Ke = Stiffness_matrix_M(E,nu,k,lex,ley,h);
Fe = Thermal_stress_M(E,nu,alp,lex,ley,h,dT);

ntot = (2*N1+1)*(N2+1)+(N1+1)*N2;

F = zeros(ntot*3,1);

II = zeros(N1*N2*24^2,1);
JJ = zeros(N1*N2*24^2,1);
SSK = zeros(N1*N2*24^2,1);
index=0;
for J = 1:N2
    for I = 1:N1
        DofLEe = 6*(I-1)+3*(3*N1+2)*(J-1)+ ...
                [(1:3) (7:9) 3*(3*N1+2)+(7:9) 3*(3*N1+2)+(1:3) ...
                 (4:6) 3*(2*N1+1)-3*(I-1)+(4:6) 3*(3*N1+2)+(4:6) 3*(2*N1+1)-3*(I-1)+(1:3)];
        F(DofLEe,1) = F(DofLEe,1) +Fe;
        for tII = 1:24
            for tJJ = 1:24
                index=index+1;
                II(index,1) = DofLEe(1,tII);
                JJ(index,1) = DofLEe(1,tJJ);
                SSK(index) = Ke(tII,tJJ);
            end
        end
    end
end

K = sparse(II,JJ,SSK,ntot*3,ntot*3);

%% kinematic constraints
DOFb=zeros(2*(N1+N2+1)+1,1);
select_w=zeros(ntot,1);

DOFb(1:3)=1:3;
select_w(1:2*N1+1)=[1 1+(1:2:(4*N1-1))];
for I = 1:2*N1+1
    DOFb(I+3) = 3*I+2;
end

for J = 1:N2
    DOFb(2*N1+3+2*(J-1)+(1:2)) = 3*(3*N1+2)*(J-1)+3*(2*N1+1)+3*[1 N1+2];
    select_w(2*N1+1+(2+3*N1)*(J-1)+(1:(2+3*N1))) = ...
       (4+9*N1)*(J-1)+(4*N1+1)+[1 (3:3:3*N1+2) 3*N1+3 3*N1+2+2+(1:3:6*N1)];
end

DOFg=1:(3*ntot);
DOFI=setdiff(DOFg,DOFb);

Ku=K(DOFI,DOFI);
Fu=F(DOFI,1);

u=Ku\Fu;
u_g=[0;u];
w=u_g(select_w);
[~,N,~] = Shape_function_M(0,0);
Nr=[N(1) N(5) N(2) N(8) N(6) N(4) N(7) N(3)];
a=0;
u_x=zeros(N1*N2,1);
for J=1:N2
    for I=1:N1
        a=a+1;
        Dof = (3*N1+2)*(J-1)+ ...
              [2*(I-1)+(1:3) (I-1)+(2*N1+1)+(1:2) 2*(I-1)+(3*N1+2)+(1:3)];
        u_x(a)=Nr*w(Dof); 
    end
end

w_plot=zeros((2*N1+1)*(2*N2+1),1);
a=0;
Do=zeros(N1*N2,1);
for J=1:N2
    for I=1:N1
        a=a+1;
        Do(a) = 2*(I-1)+(2*N1+1)*2*(J-1)+(2*N1+3);
    end
end
w_plot(Do,1)=u_x;
Di=setdiff(1:(2*N1+1)*(2*N2+1),Do);
w_plot(Di,1)=w;

X=0:lex/2:Lx;
Y=0:ley/2:Ly;

[XX,YY]=meshgrid(X,Y);

ZZ=zeros(size(XX));

figure(1);

q_plot=zeros(2*N1+1,2*N2+1);
q_plot(1:(2*N1+1),1:(2*N2+1))=reshape(-w_plot,2*N1+1,2*N2+1);

surf(XX,YY,ZZ)
hold on
surf(XX,YY,q_plot)
axis([0 Lx 0 Ly])
title("Displacement |Uz|(m) - Mindlin plate")
colorbar
shading interp

u_M=zeros((N1+1)*(N2+1),1);

for J = 1:N2+1
    for I = 1:N1+1
        u_M(I+(J-1)*(N1+1)) = w_plot(2*(I-1)+1+2*(J-1)*(2*N1+1));
    end
end
