function[M] = Mass_matrix_M(rho,lex,ley,h)

nodeCoor = [0 0; lex 0; lex ley; 0 ley; lex/2 0; lex ley/2; lex/2 ley; 0 ley/2];

Ie = [1 0 0; 0 h^2/12 0; 0 0 h^2/12];

xi = [-1/sqrt(3) 1/sqrt(3)];
eta = xi;
c1 = [1 1];
c2 = c1;
M = zeros(24,24);
for I=1:2
    for J = 1:2
        [N,~,~] = Shape_function_M(xi(I),eta(J));
        [Jac,~] = Jacobian_M(nodeCoor, xi(I), eta(J));
        M = c1(I)*c2(J)*N.'*Ie*N*det(Jac) + M;
    end
end

M = rho*h*M;

end
