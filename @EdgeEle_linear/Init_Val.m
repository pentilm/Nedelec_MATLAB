function Init_Val( obj )
% deal with the initial value of all functions by using mass matrix
% find E0, J0, H0, K0 by using initial value condition
% use the function obj.Ex_exa, obj.Ey_exa, obj.Jx_exa, obj.Jy_exa,
% obj.H_exa, obj.K_exa at t=0 as the I.V.

raw_E = zeros(obj.nedg,1);
raw_J = zeros(obj.nedg,1);
raw_H = zeros(obj.nele,1);
raw_K = zeros(obj.nele,1);

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
    for j = 1:obj.side_ele
        ed1 = obj.el2ed(j,i);
        % for E, J
        switch j
        % (f,\phi)
            case {1,3}
                fcn_E = @(x,y)obj.Ex_exa(x,y,obj.at).*Basis_fcn1(x,y,i,j,obj);
                fcn_J = @(x,y)obj.Jx_exa(x,y,obj.at).*Basis_fcn1(x,y,i,j,obj);
            case {2,4}
                fcn_E = @(x,y)obj.Ey_exa(x,y,obj.at).*Basis_fcn1(x,y,i,j,obj);
                fcn_J = @(x,y)obj.Jy_exa(x,y,obj.at).*Basis_fcn1(x,y,i,j,obj);
        end
        int_val_E = Quad_Ele_Rect(@(x,y)fcn_E(x,y),[xae,yae],[xbe,ybe]);
        int_val_J = Quad_Ele_Rect(@(x,y)fcn_J(x,y),[xae,yae],[xbe,ybe]);
        raw_E(ed1) = raw_E(ed1) + int_val_E*obj.edori(i,j); % be careful for the direction of edges
        raw_J(ed1) = raw_J(ed1) + int_val_J*obj.edori(i,j);
    end
	% for H, K
    fcn_H = @(x,y)obj.H_exa(x,y,obj.at);
    fcn_K = @(x,y)obj.K_exa(x,y,obj.at);
    raw_H(i) = Quad_Ele_Rect(@(x,y)fcn_H(x,y),[xae,yae],[xbe,ybe]);
    raw_K(i) = Quad_Ele_Rect(@(x,y)fcn_K(x,y),[xae,yae],[xbe,ybe]);
% % another chioce: use mid point as the intial value
% raw_H(i) = fcn_H(0.5*(xae+xbe),0.5*(yae+ybe));
% raw_K(i) = fcn_K(0.5*(xae+xbe),0.5*(yae+ybe));
end
obj.E0 = obj.matM\raw_E(obj.eint);
obj.J0 = obj.matM\raw_J(obj.eint);
obj.H0 = raw_H./obj.ele_area;
obj.K0 = raw_K./obj.ele_area;
% % another chioce: use mid point as the intial value
% obj.H0 = raw_H;
% obj.K0 = raw_K;
end

