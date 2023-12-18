function[Bf, Bs] = Matrix_der_M(nodeCoor, xi,eta)

[~,N_row,~] = Shape_function_M(xi,eta);

[~,XYder] = Jacobian_M(nodeCoor,xi,eta);

Bf(1,3:3:24) = XYder(1,:);
Bf(2,2:3:23) = -XYder(2,:);
Bf(3,2:3:23) = -XYder(1,:);
Bf(3,3:3:24) = XYder(2,:);

Bs(1,1:3:22) = XYder(1,:);
Bs(1,3:3:24) = N_row;
Bs(2,1:3:22) = XYder(2,:);
Bs(2,2:3:23) = -N_row;

end
