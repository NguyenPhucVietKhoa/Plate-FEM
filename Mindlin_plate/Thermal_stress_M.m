function[F] = Thermal_stress_M(E,nu,alp,lex,ley,h,dT)

nodeCoor = [0 0; lex 0; lex ley; 0 ley; lex/2 0; lex ley/2; lex/2 ley; 0 ley/2];

Hf = E*h^2/(12*(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

eps_ther = alp*dT*[1; 1; 0];

xi = [-1/sqrt(3) 1/sqrt(3)];
eta = xi;

F = zeros(24,1);

for I=1:2
    for J = 1:2
        [Bf,~] = Matrix_der_M(nodeCoor,xi(I),eta(J));
        [Jac,~] = Jacobian_M(nodeCoor, xi(I), eta(J));
        F = Bf.'*Hf*eps_ther*det(Jac) + F;
    end
end

end