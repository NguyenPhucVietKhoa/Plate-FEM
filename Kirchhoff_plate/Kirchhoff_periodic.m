%%
clear
close all
%clc

%Material parameters
%substrat
E_s=1.3e11;
nu_s=0.27;
alp_s=2.8e-6;

%Couches
E_c=69e9;
nu_c=0.33;
alp_c=24e-6;
%Chargement
dT = 100;
%Maillage
%nombre des sous-strucutres
Nx = 5;
Ny = 5;
%nombre des éléments par sous-structures
Ns1 = 10;
Ns2 = 10;
ns1 = Ns1+1;
ns2 = Ns2+1;
%Nombre des éléments totals
N1 = Nx*Ns1;
N2 = Ny*Ns2;

%Dimensions
Lx=1;
Ly=1;
h_s=0.001; % Substrat
h_c=0.003; % Couches

lex = Lx/N1;
ley = Ly/N2;

%Nombre des éléments couches dans chaque sous-structure
Nc1=4;
Nc2=4;
NN1 = (Ns1-Nc1)/2; 
NN2 = (Ns2-Nc2)/2;

Ke_s = Stiffness_matrix_K(E_s,nu_s,lex,ley,h_s);
Fe_s = Thermal_stress_K(E_s,nu_s,alp_s,lex,ley,h_s,dT);

Ke_c = Stiffness_matrix_K(E_c,nu_c,lex,ley,h_c);
Fe_c = Thermal_stress_K(E_c,nu_c,alp_c,lex,ley,h_c,dT);
%number of nodes along x1:
n1=N1+1;
%number of nodes along x2:
n2=N2+1;

II = zeros(N1*N2*12^2,1);
JJ = zeros(N1*N2*12^2,1);
SSK = zeros(N1*N2*12^2,1);
index=0;

F=zeros(n1*n2*3,1);
for k=1:Ny
    for o=1:Nx
        for J = 1:Ns2
            for I = 1:Ns1
                DofLEe = 3*Ns1*(o-1)+3*n1*(ns2-1)*(k-1)+ ...
                    [3*(J-1)*n1+3*(I-1)+(1:6) 3*J*n1+3*(I-1)+(4:6) 3*J*n1+3*(I-1)+(1:3)];

                if I > NN1 && J > NN2 && I <= Ns1-NN1 && J <= Ns2-NN1
                    F(DofLEe,1) = F(DofLEe,1) + Fe_c;
                    for tII = 1:12
                        for tJJ = 1:12
                            index=index+1;
                            II(index,1) = DofLEe(1,tII);
                            JJ(index,1) = DofLEe(1,tJJ);
                            SSK(index) = Ke_c(tII,tJJ);
                        end
                    end
                else
                    F(DofLEe,1) = F(DofLEe,1) + Fe_s;
                    for tII = 1:12
                        for tJJ = 1:12
                            index=index+1;
                            II(index,1) = DofLEe(1,tII);
                            JJ(index,1) = DofLEe(1,tJJ);
                            SSK(index) = Ke_s(tII,tJJ);
                        end
                    end
                end
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

u_glo = zeros(3*n1*n2,1);
u_glo(DOFI,1) = u;

DOF_coins = zeros(3*(Nx+1)*(Ny+1),1);
for J=1:Ny+1
    for I = 1:Nx+1
        DOF_coins(3*(I-1)+3*(Nx+1)*(J-1)+(1:3)) = 3*n1*(ns2-1)*(J-1)+3*(ns1-1)*(I-1)+(1:3);
    end
end
u_classique = u_glo(DOF_coins,1);

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

shading interp