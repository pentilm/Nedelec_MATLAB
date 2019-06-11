function Marching( obj )
cst1 = 1+obj.dt^2*obj.we^2/(4+2*obj.dt*obj.ge);
cst2 = 1+obj.dt^2*obj.wm^2/(4+2*obj.dt*obj.gm);
mat_A = cst1*obj.matM;
mat_B = 0.5*obj.dt*obj.matBM;
% mat_D = cst2*obj.matC;
mat_D_inv = 1/cst2*obj.matC_inv;
mat_E = mat_A + mat_B*mat_D_inv*mat_B';
tic;
for it = 2:obj.Nt
ct = obj.dt*(it-1);
obj.Form_RHS;
Form_RHS_tmp( obj,ct );
vec_f = obj.RHS + obj.RHS_tmp;  % whole RHS vector
vec_f1 = vec_f(1:obj.ienum);
vec_f2 = vec_f(obj.ienum+1:end);
obj.E1 = mat_E\(vec_f1 + mat_B*mat_D_inv*vec_f2);
obj.H1 = mat_D_inv*(vec_f2-mat_B'*obj.E1);
UpdateJK( obj );
obj.E0 = obj.E1;
obj.H0 = obj.H1;
obj.J0 = obj.J1;
obj.K0 = obj.K1;
end
toc;

end

