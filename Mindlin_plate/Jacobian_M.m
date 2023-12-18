function [Jacobian, XYder] = Jacobian_M(nodeCoor,xi,eta)

[~,~,derN] = Shape_function_M(xi,eta);

Jacobian(1,1) = derN(1,:)*nodeCoor(:,1);
Jacobian(1,2) = derN(2,:)*nodeCoor(:,1);
Jacobian(2,1) = derN(1,:)*nodeCoor(:,2);
Jacobian(2,2) = derN(2,:)*nodeCoor(:,2);

XYder = inv(Jacobian).'*derN;

end