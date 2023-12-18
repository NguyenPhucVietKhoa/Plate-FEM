function [Jacobian] = Jacobian_K(xi,eta,lex,ley)

nodeCoor = [0 0; lex 0; lex ley; 0 ley];

%Q4=1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];

derQ4 = 1/4*[-(1-eta) 1-eta 1+eta -(1+eta); ...
             -(1-xi) -(1+xi) 1+xi 1-xi];

Jacobian(1,1) = derQ4(1,:)*nodeCoor(:,1);
Jacobian(1,2) = derQ4(2,:)*nodeCoor(:,1);
Jacobian(2,1) = derQ4(1,:)*nodeCoor(:,2);
Jacobian(2,2) = derQ4(2,:)*nodeCoor(:,2);

end