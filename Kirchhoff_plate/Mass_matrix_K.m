function[M] = Mass_matrix_K(rho,lex,ley,h)

xi = [-0.906 -0.538 0 0.538 0.906];
eta = xi;
c1 = [0.237 0.479 0.569 0.479 0.237];
c2 = c1;
M = zeros(12,12);

for I=1:5
    for J=1:5
        [N,~] = Shape_function_K(xi(I),eta(J),lex,ley);
        [Ja] = Jacobian_K(xi(I),eta(J),lex,ley);
        M = c1(I)*c2(J)*(N.'*N)*det(Ja) + M;
    end
end

M = rho*h*M;
end