function[F] = Thermal_stress_K(E,nu,alp,lex,ley,h,dT)

Hf = E*h^2/(12*(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

eps_ther = alp*dT*[1; 1; 0];

xi = [-1/sqrt(3) 1/sqrt(3)];
eta = xi;

F = zeros(12,1);

for I=1:2
    for J=1:2
        [B] = Matrix_der_K(xi(I),eta(J),lex,ley);
        [Ja] = Jacobian_K(xi(I),eta(J),lex,ley);
        F = B.'*Hf*eps_ther*det(Ja) + F;
    end
end

end
