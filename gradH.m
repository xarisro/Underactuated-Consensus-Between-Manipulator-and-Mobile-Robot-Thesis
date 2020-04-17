%gradH.m
%Calculating the gradient of H, dH/dq
function grad = gradH(dofs,J,robot,q,end_eff_name)
dq = 0.002;
ndof = length(dofs);
Jtemp = zeros(6,ndof);
grad = zeros(ndof,1);

initH = Hfunc(J);
for i=1:ndof
    qtemp = q;
    qtemp(dofs(i)) = qtemp(dofs(i)) + dq;
    
    fullJac = geometricJacobian(robot,qtemp,end_eff_name);
    for j=1:ndof
        Jtemp(:,j) = fullJac(:,dofs(j));
    end
    
    grad(i) = (Hfunc(Jtemp) - initH)/dq; %Arithmetical derivation
end

end

