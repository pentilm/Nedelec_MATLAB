function Post_Process( obj )
% get the function values at the center of every elements
% compute error of L^\infty and L^2
% for rectangle partition, just use matrix
% for triangle partition, use column vector and then assemble
% err_inf, err_l2: Ex, Ey, Jx, Jy, H, K

% 1-D vector for numerical solution. Used for all partition cases
Ex_vec = zeros(obj.nele,1);
Ey_vec = zeros(obj.nele,1);
Jx_vec = zeros(obj.nele,1);
Jy_vec = zeros(obj.nele,1);
H_vec = zeros(obj.nele,1);
K_vec = zeros(obj.nele,1);
% 1-D vector for exact solution. Used for all partition cases
Ex_ex_vec = zeros(obj.nele,1);
Ey_ex_vec = zeros(obj.nele,1);
Jx_ex_vec = zeros(obj.nele,1);
Jy_ex_vec = zeros(obj.nele,1);
H_ex_vec = zeros(obj.nele,1);
K_ex_vec = zeros(obj.nele,1);

% intialize numerical solution matrix on rectangle partition
obj.Ex_num = zeros(obj.Ny-1,obj.Nx-1);
obj.Ey_num = zeros(obj.Ny-1,obj.Nx-1);
obj.Jx_num = zeros(obj.Ny-1,obj.Nx-1);
obj.Jy_num = zeros(obj.Ny-1,obj.Nx-1);
obj.H_num = zeros(obj.Ny-1,obj.Nx-1);
obj.K_num = zeros(obj.Ny-1,obj.Nx-1);

%% expand the E1 and J1 to the coefficients for whole edges with boundary ones
% for boundary edges, set coefficients as 0
solvec_E = zeros(obj.nedg,1);
solvec_E(obj.eint) = obj.E1;
solvec_J = zeros(obj.nedg,1);
solvec_J(obj.eint) = obj.J1;

for i = 1:obj.nele   % loop for all elements
% coordinates of this element from 1st node & 3rd node
%
%    (xae,ybe)----(xbe,ybe)
%       |           |
%       |           |
%    (xae,yae)----(xbe,yae)
%
   xae = obj.no2xy(1,obj.el2no(1,i));  xbe = obj.no2xy(1,obj.el2no(3,i));
   yae = obj.no2xy(2,obj.el2no(1,i));  ybe = obj.no2xy(2,obj.el2no(3,i));
   mid_x = 0.5*(xae+xbe);  mid_y = 0.5*(yae+ybe);
%% write exact and numerical solution at mid point of each element
%% Ex
Ex_vec(i) = obj.edori(i,1)*solvec_E(obj.el2ed(1,i))*Basis_fcn1(mid_x,mid_y,i,1,obj) ...
          + obj.edori(i,3)*solvec_E(obj.el2ed(3,i))*Basis_fcn1(mid_x,mid_y,i,3,obj);
Ex_ex_vec(i) = obj.Ex_exa(mid_x,mid_y,obj.bt);
%% Ey
Ey_vec(i) = obj.edori(i,2)*solvec_E(obj.el2ed(2,i))*Basis_fcn1(mid_x,mid_y,i,2,obj) ...
           +obj.edori(i,4)*solvec_E(obj.el2ed(4,i))*Basis_fcn1(mid_x,mid_y,i,4,obj);
Ey_ex_vec(i) = obj.Ey_exa(mid_x,mid_y,obj.bt);
%% Jx
Jx_vec(i) = obj.edori(i,1)*solvec_J(obj.el2ed(1,i))*Basis_fcn1(mid_x,mid_y,i,1,obj) ...
          + obj.edori(i,3)*solvec_J(obj.el2ed(3,i))*Basis_fcn1(mid_x,mid_y,i,3,obj);
Jx_ex_vec(i) = obj.Jx_exa(mid_x,mid_y,obj.bt);
%% Jy
Jy_vec(i) = obj.edori(i,2)*solvec_J(obj.el2ed(2,i))*Basis_fcn1(mid_x,mid_y,i,2,obj) ...
           +obj.edori(i,4)*solvec_J(obj.el2ed(4,i))*Basis_fcn1(mid_x,mid_y,i,4,obj);
Jy_ex_vec(i) = obj.Jy_exa(mid_x,mid_y,obj.bt);
%% H
H_vec(i) = obj.H1(i);
H_ex_vec(i) = obj.H_exa(mid_x,mid_y,obj.bt);
%% K
K_vec(i) = obj.K1(i);
K_ex_vec(i) = obj.K_exa(mid_x,mid_y,obj.bt);

%% L^2 error
%% Ex
fcn_tmp_Ex = @(x,y)(obj.edori(i,1)*solvec_E(obj.el2ed(1,i))*Basis_fcn1(x,y,i,1,obj) ...
          + obj.edori(i,3)*solvec_E(obj.el2ed(3,i))*Basis_fcn1(x,y,i,3,obj) - obj.Ex_exa(x,y,obj.bt)).^2;
int_val_Ex = Quad_Ele_Rect(@(x,y)fcn_tmp_Ex(x,y),[xae,yae],[xbe,ybe]);
obj.err_l2(1) = obj.err_l2(1) + int_val_Ex;
%% Ey
fcn_tmp_Ey = @(x,y)(obj.edori(i,2)*solvec_E(obj.el2ed(2,i))*Basis_fcn1(x,y,i,2,obj) ...
           +obj.edori(i,4)*solvec_E(obj.el2ed(4,i))*Basis_fcn1(x,y,i,4,obj) - obj.Ey_exa(x,y,obj.bt)).^2;
int_val_Ey = Quad_Ele_Rect(@(x,y)fcn_tmp_Ey(x,y),[xae,yae],[xbe,ybe]);
obj.err_l2(2) = obj.err_l2(2) + int_val_Ey;
%% Jx
fcn_tmp_Jx = @(x,y)(obj.edori(i,1)*solvec_J(obj.el2ed(1,i))*Basis_fcn1(x,y,i,1,obj) ...
          + obj.edori(i,3)*solvec_J(obj.el2ed(3,i))*Basis_fcn1(x,y,i,3,obj) - obj.Jx_exa(x,y,obj.bt)).^2;
int_val_Jx = Quad_Ele_Rect(@(x,y)fcn_tmp_Jx(x,y),[xae,yae],[xbe,ybe]);
obj.err_l2(3) = obj.err_l2(3) + int_val_Jx;
%% Jy
fcn_tmp_Jy = @(x,y)(obj.edori(i,2)*solvec_J(obj.el2ed(2,i))*Basis_fcn1(x,y,i,2,obj) ...
           +obj.edori(i,4)*solvec_J(obj.el2ed(4,i))*Basis_fcn1(x,y,i,4,obj) - obj.Jy_exa(x,y,obj.bt)).^2;
int_val_Jy = Quad_Ele_Rect(@(x,y)fcn_tmp_Jy(x,y),[xae,yae],[xbe,ybe]);
obj.err_l2(4) = obj.err_l2(4) + int_val_Jy;
%% H
fcn_tmp_H = @(x,y)(obj.H1(i)-obj.H_exa(x,y,obj.bt)).^2;
int_val_H = Quad_Ele_Rect(@(x,y)fcn_tmp_H(x,y),[xae,yae],[xbe,ybe]);
obj.err_l2(5) = obj.err_l2(5) + int_val_H;
%% K
fcn_tmp_K = @(x,y)(obj.K1(i)-obj.K_exa(x,y,obj.bt)).^2;
int_val_K = Quad_Ele_Rect(@(x,y)fcn_tmp_K(x,y),[xae,yae],[xbe,ybe]);
obj.err_l2(6) = obj.err_l2(6) + int_val_K;   
end
%% deal with L^2 error
obj.err_l2 = sqrt(obj.err_l2);
%% deal with L^infty error
obj.err_inf(1) = max(abs(Ex_vec-Ex_ex_vec));
obj.err_inf(2) = max(abs(Ey_vec-Ey_ex_vec));
obj.err_inf(3) = max(abs(Jx_vec-Jx_ex_vec));
obj.err_inf(4) = max(abs(Jy_vec-Jy_ex_vec));
obj.err_inf(5) = max(abs(H_vec-H_ex_vec));
obj.err_inf(6) = max(abs(K_vec-K_ex_vec));
%% convert 1D vector to 2D matrix
nelex = obj.Nx-1;
for j = 1:(obj.Ny-1)
    for i = 1:(obj.Nx-1)
        obj.Ex_num(j,i) = Ex_vec(nelex*(j-1)+i);
        obj.Ey_num(j,i) = Ey_vec(nelex*(j-1)+i);
        obj.Jx_num(j,i) = Jx_vec(nelex*(j-1)+i);
        obj.Jy_num(j,i) = Jy_vec(nelex*(j-1)+i);
        obj.H_num(j,i) = H_vec(nelex*(j-1)+i);
        obj.K_num(j,i) = K_vec(nelex*(j-1)+i);
    end
end

end

