function[B] = Matrix_der_K(xi,eta,lex,ley)

[~,derN] = Shape_function_K(xi,eta,lex,ley);

[Jac] = Jacobian_K(xi, eta, lex,ley);

J_invT = inv(Jac).';

B(1,:) = -J_invT(1,1)^2*derN(3,:);
B(2,:) = -J_invT(2,2)^2*derN(4,:);
B(3,:) = -2*J_invT(1,1)*J_invT(2,2)*derN(5,:);

end