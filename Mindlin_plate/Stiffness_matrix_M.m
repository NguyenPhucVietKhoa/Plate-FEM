function[K] = Stiffness_matrix_M(E,nu,k,lex,ley,h)

nodeCoor = [0 0; lex 0; lex ley; 0 ley; lex/2 0; lex ley/2; lex/2 ley; 0 ley/2];

Hf = E*h^2/(12*(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

Hs = E*k/(2*(1+nu))*eye(2,2);

xi = [-1/sqrt(3) 1/sqrt(3)];
eta = xi;

Kf = zeros(24,24);
Ks = zeros(24,24);
for I=1:2
    for J = 1:2
        [Bf,Bs] = Matrix_der_M(nodeCoor,xi(I),eta(J));
        [Jac,~] = Jacobian_M(nodeCoor, xi(I), eta(J));
        Kf = Bf.'*Hf*Bf*det(Jac) + Kf;
        Ks = Bs.'*Hs*Bs*det(Jac) + Ks;
    end
end

K =h*(Kf + Ks);

end
