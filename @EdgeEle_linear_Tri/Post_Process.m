function Post_Process( obj )
% get the function values at the center of every elements
% compute error of L^\infty and L^2
% for rectangle partition, just use matrix
% for triangle partition, use column vector and then assemble
% err_inf, err_l2: Ex, Ey, Jx, Jy, H, K

% 1-D vector for numerical solution. Used for all partition cases
Ex_vec = zeros(obj.npt,1);
Ey_vec = zeros(obj.npt,1);
Jx_vec = zeros(obj.npt,1);
Jy_vec = zeros(obj.npt,1);
H_vec = zeros(obj.nele,1);
K_vec = zeros(obj.nele,1);
% 1-D vector for exact solution. Used for all partition cases
Ex_ex_vec = zeros(obj.npt,1);
Ey_ex_vec = zeros(obj.npt,1);
Jx_ex_vec = zeros(obj.npt,1);
Jy_ex_vec = zeros(obj.npt,1);
H_ex_vec = zeros(obj.nele,1);
K_ex_vec = zeros(obj.nele,1);
% auxiliary arrays used to count how many times to visit XX_vec vector
Ex_count = zeros(obj.npt,1);
Ey_count = zeros(obj.npt,1);
Jx_count = zeros(obj.npt,1);
Jy_count = zeros(obj.npt,1);

%% expand the E1 and J1 to the coefficients for whole edges with boundary ones
% for boundary edges, set coefficients as 0
solvec_E = zeros(obj.nedg,1);
solvec_E(obj.eint) = obj.E1;
solvec_J = zeros(obj.nedg,1);
solvec_J(obj.eint) = obj.J1;

for i = 1:obj.nele   % loop for all elements
% coordinates of this element from 1st node & 3rd node
%             (x1,y1)
%               / \
%              /   \
%             /     \
%        (x2,y2)---(x3,y3)    
%
   x1 = obj.no2xy(1,obj.el2no(1,i));  y1 = obj.no2xy(2,obj.el2no(1,i));
   x2 = obj.no2xy(1,obj.el2no(2,i));  y2 = obj.no2xy(2,obj.el2no(2,i));
   x3 = obj.no2xy(1,obj.el2no(3,i));  y3 = obj.no2xy(2,obj.el2no(3,i));
   p = [x1 , x2 , x3; y1 , y2 , y3];
%% write exact and numerical solution at mid point of each element
%% to test the superconvergence of triangle element, find the error on every point, and then compute their average value
for j = 1:obj.side_ele
	% point of the vertex
    pt = obj.el2no(j,i);
	pt_x = obj.no2xy(1,pt);
	pt_y = obj.no2xy(2,pt);
%% Ex
Ex_vec(pt) = Ex_vec(pt) ...
          + obj.edori(i,1)*solvec_E(obj.el2ed(1,i))*Basis_fcn1(pt_x,pt_y,i,1,1,obj) ...
          + obj.edori(i,2)*solvec_E(obj.el2ed(2,i))*Basis_fcn1(pt_x,pt_y,i,2,1,obj) ...
          + obj.edori(i,3)*solvec_E(obj.el2ed(3,i))*Basis_fcn1(pt_x,pt_y,i,3,1,obj);
Ex_ex_vec(pt) = obj.Ex_exa(pt_x,pt_y,obj.bt);
Ex_count(pt) = Ex_count(pt) + 1;
%% Ey
Ey_vec(pt) = Ey_vec(pt) ...
           + obj.edori(i,1)*solvec_E(obj.el2ed(1,i))*Basis_fcn1(pt_x,pt_y,i,1,2,obj) ...
           + obj.edori(i,2)*solvec_E(obj.el2ed(2,i))*Basis_fcn1(pt_x,pt_y,i,2,2,obj) ...
           + obj.edori(i,3)*solvec_E(obj.el2ed(3,i))*Basis_fcn1(pt_x,pt_y,i,3,2,obj);
Ey_ex_vec(pt) = obj.Ey_exa(pt_x,pt_y,obj.bt);
Ey_count(pt) = Ey_count(pt) + 1;
%% Jx
Jx_vec(pt) = Jx_vec(pt) ...
          + obj.edori(i,1)*solvec_J(obj.el2ed(1,i))*Basis_fcn1(pt_x,pt_y,i,1,1,obj) ...
          + obj.edori(i,2)*solvec_J(obj.el2ed(2,i))*Basis_fcn1(pt_x,pt_y,i,2,1,obj) ...
          + obj.edori(i,3)*solvec_J(obj.el2ed(3,i))*Basis_fcn1(pt_x,pt_y,i,3,1,obj);
Jx_ex_vec(pt) = obj.Jx_exa(pt_x,pt_y,obj.bt);
Jx_count(pt) = Jx_count(pt) + 1;
%% Jy
Jy_vec(pt) = Jy_vec(pt) ...
          + obj.edori(i,1)*solvec_J(obj.el2ed(1,i))*Basis_fcn1(pt_x,pt_y,i,1,2,obj) ...
          + obj.edori(i,2)*solvec_J(obj.el2ed(2,i))*Basis_fcn1(pt_x,pt_y,i,2,2,obj) ...
          + obj.edori(i,3)*solvec_J(obj.el2ed(3,i))*Basis_fcn1(pt_x,pt_y,i,3,2,obj);
Jy_ex_vec(pt) = obj.Jy_exa(pt_x,pt_y,obj.bt);
Jy_count(pt) = Jy_count(pt) + 1;
end
%% H
H_vec(i) = obj.H1(i);
H_ex_vec(i) = obj.H_exa(obj.px_mid(i),obj.py_mid(i),obj.bt);
%% K
K_vec(i) = obj.K1(i);
K_ex_vec(i) = obj.K_exa(obj.px_mid(i),obj.py_mid(i),obj.bt);
%% L^2 error
%% Ex
fcn_tmp_Ex = @(x,y)(obj.edori(i,1)*solvec_E(obj.el2ed(1,i))*Basis_fcn1(x,y,i,1,1,obj) ...
          + obj.edori(i,2)*solvec_E(obj.el2ed(2,i))*Basis_fcn1(x,y,i,2,1,obj) ...
          + obj.edori(i,3)*solvec_E(obj.el2ed(3,i))*Basis_fcn1(x,y,i,3,1,obj) ...
          - obj.Ex_exa(x,y,obj.bt)).^2;
int_val_Ex = Quad_Ele_Tri(@(x,y)fcn_tmp_Ex(x,y),p);
obj.err_l2(1) = obj.err_l2(1) + int_val_Ex;
%% Ey
fcn_tmp_Ey = @(x,y)(obj.edori(i,1)*solvec_E(obj.el2ed(1,i))*Basis_fcn1(x,y,i,1,2,obj) ...
           +obj.edori(i,2)*solvec_E(obj.el2ed(2,i))*Basis_fcn1(x,y,i,2,2,obj) ...
           +obj.edori(i,3)*solvec_E(obj.el2ed(3,i))*Basis_fcn1(x,y,i,3,2,obj) ...
           - obj.Ey_exa(x,y,obj.bt)).^2;
int_val_Ey = Quad_Ele_Tri(@(x,y)fcn_tmp_Ey(x,y),p);
obj.err_l2(2) = obj.err_l2(2) + int_val_Ey;
%% Jx
fcn_tmp_Jx = @(x,y)(obj.edori(i,1)*solvec_J(obj.el2ed(1,i))*Basis_fcn1(x,y,i,1,1,obj) ...
          + obj.edori(i,2)*solvec_J(obj.el2ed(2,i))*Basis_fcn1(x,y,i,2,1,obj) ...
          + obj.edori(i,3)*solvec_J(obj.el2ed(3,i))*Basis_fcn1(x,y,i,3,1,obj) ...
          - obj.Jx_exa(x,y,obj.bt)).^2;
int_val_Jx = Quad_Ele_Tri(@(x,y)fcn_tmp_Jx(x,y),p);
obj.err_l2(3) = obj.err_l2(3) + int_val_Jx;
%% Jy
fcn_tmp_Jy = @(x,y)(obj.edori(i,1)*solvec_J(obj.el2ed(1,i))*Basis_fcn1(x,y,i,1,2,obj) ...
          + obj.edori(i,2)*solvec_J(obj.el2ed(2,i))*Basis_fcn1(x,y,i,2,2,obj) ...
          + obj.edori(i,3)*solvec_J(obj.el2ed(3,i))*Basis_fcn1(x,y,i,3,2,obj) ...
          - obj.Jy_exa(x,y,obj.bt)).^2;
int_val_Jy = Quad_Ele_Tri(@(x,y)fcn_tmp_Jy(x,y),p);
obj.err_l2(4) = obj.err_l2(4) + int_val_Jy;
%% H
fcn_tmp_H = @(x,y)(obj.H1(i)-obj.H_exa(x,y,obj.bt)).^2;
int_val_H = Quad_Ele_Tri(@(x,y)fcn_tmp_H(x,y),p);
obj.err_l2(5) = obj.err_l2(5) + int_val_H;
%% K
fcn_tmp_K = @(x,y)(obj.K1(i)-obj.K_exa(x,y,obj.bt)).^2;
int_val_K = Quad_Ele_Tri(@(x,y)fcn_tmp_K(x,y),p);
obj.err_l2(6) = obj.err_l2(6) + int_val_K;   
end
%% calculate the average on each edge
Ex_vec = Ex_vec./Ex_count;
Ey_vec = Ey_vec./Ey_count;
Jx_vec = Jx_vec./Jx_count;
Jy_vec = Jy_vec./Jy_count;
%% deal with L^2 error
obj.err_l2 = sqrt(obj.err_l2);
%% deal with L^infty error
% mark=ones(obj.nedg,1);
nn=1;
for i =1:obj.npt
    p0 = obj.no2xy(:,i);
    p = [0.5;0.5];
    d = sqrt(sum((p0-p).^2));
    if(d>0.1)
        continue;
    else
        mark(nn) = i;
        nn=nn+1;
    end
end
% obj.err_inf(1) = max(abs(Ex_vec(obj.eint)-Ex_ex_vec(obj.eint)));
% obj.err_inf(2) = max(abs(Ey_vec(obj.eint)-Ey_ex_vec(obj.eint)));
% obj.err_inf(3) = max(abs(Jx_vec(obj.eint)-Jx_ex_vec(obj.eint)));
% obj.err_inf(4) = max(abs(Jy_vec(obj.eint)-Jy_ex_vec(obj.eint)));
obj.err_inf(1) = max(abs(Ex_vec(mark)-Ex_ex_vec(mark)));
obj.err_inf(2) = max(abs(Ey_vec(mark)-Ey_ex_vec(mark)));
obj.err_inf(3) = max(abs(Jx_vec(mark)-Jx_ex_vec(mark)));
obj.err_inf(4) = max(abs(Jy_vec(mark)-Jy_ex_vec(mark)));

obj.err_inf(5) = max(abs(H_vec-H_ex_vec));
obj.err_inf(6) = max(abs(K_vec-K_ex_vec));
%% for triangle mesh, just get the 1-D solution vector
obj.Ex_num = Ex_vec;
obj.Ey_num = Ey_vec;
obj.Jx_num = Jx_vec;
obj.Jy_num = Jy_vec;
obj.H_num = H_vec;
obj.K_num = K_vec;

end

