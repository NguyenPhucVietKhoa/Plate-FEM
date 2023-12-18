function[K] = Stiffness_matrix_K(E,nu,lex,ley,h)

Hf = E*h^2/(12*(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

xi = [-sqrt(3/5) 0 sqrt(3/5)];
eta = xi;
c1 = [5/9 8/9 5/9];
c2 = c1;

K = zeros(12,12);

for I=1:3
    for J=1:3
        [B] = Matrix_der_K(xi(I),eta(J),lex,ley);
        [Ja] = Jacobian_K(xi(I),eta(J),lex,ley);
        K = h*c1(I)*c2(J)*B.'*Hf*B*det(Ja) + K;
    end
end

end