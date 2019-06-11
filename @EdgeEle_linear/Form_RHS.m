function Form_RHS( obj )
% assemble the time independent RHS vector, run for every time loop
obj.RHS = zeros(obj.num_unk,1);
% for E
% update obj.RHS(1:obj.nedg)
cst1 = 1 - obj.dt^2*obj.we^2/(4+2*obj.dt*obj.ge);
cst2 = -2*obj.dt/(2+obj.dt*obj.ge);
obj.RHS(1:obj.ienum) = cst1*obj.matM*obj.E0 + 0.5*obj.dt*obj.matBM*obj.H0 + cst2*obj.matM*obj.J0;
% for H
% obj.RHS(obj.nedg+1:end)
cst1 = 1 - obj.dt^2*obj.wm^2/(4+2*obj.dt*obj.gm);
cst2 = -2*obj.dt/(2+obj.dt*obj.gm);
obj.RHS(obj.ienum+1:end) = cst1*obj.matC*obj.H0 - 0.5*obj.dt*obj.matBM'*obj.E0 + cst2*obj.matC*obj.K0; 


end

