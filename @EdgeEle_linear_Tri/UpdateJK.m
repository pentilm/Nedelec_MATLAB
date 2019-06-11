function UpdateJK( obj )
% update J and K (without boundary points)
% Since no spatial differential for J and K, just update their coefficients
% For RHS, find their projections by using vec_g3 and vec_g4

% For J
cst1 = (2-obj.dt*obj.ge)./(2+obj.dt*obj.ge);
cst2 = obj.dt*obj.we.^2./(2+obj.dt*obj.ge);
cst3 = 2*obj.dt./(2+obj.dt*obj.ge);
obj.J1 = cst1*obj.J0 + cst2*(obj.E0+obj.E1) + cst3*(obj.matM\obj.vec_g3);

% For K
cst1 = (2-obj.dt*obj.gm)./(2+obj.dt*obj.gm);
cst2 = obj.dt*obj.gm.^2./(2+obj.dt*obj.gm);
cst3 = 2*obj.dt./(2+obj.dt*obj.gm);
obj.K1 = cst1*obj.K0 + cst2*(obj.H0+obj.H1) + cst3*(obj.matC_inv*obj.vec_g4);

end

