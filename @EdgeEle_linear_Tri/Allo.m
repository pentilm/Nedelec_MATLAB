function Allo( obj )
% initialize all the arries
% computational parameters
    obj.no2xy = zeros(2,obj.npt);
	obj.el2no = zeros(obj.side_ele,obj.nele); 
	obj.edori = ones(obj.nele,obj.side_ele);
	obj.ed_id = zeros(obj.nedg,1);
	obj.ele_area = zeros(obj.nele,1);
    obj.px_mid = zeros(obj.nele,1);
    obj.py_mid = zeros(obj.nele,1);
	obj.matM = sparse(obj.nedg,obj.nedg);
	obj.matBM = sparse(obj.nedg,obj.nele);
	obj.E0 = zeros(obj.nedg,1);
	obj.H0 = zeros(obj.nele,1);
	obj.E1 = zeros(obj.nedg,1);
	obj.H1 = zeros(obj.nele,1);
	obj.J0 = zeros(obj.nedg,1);
	obj.K0 = zeros(obj.nele,1);
	obj.J1 = zeros(obj.nedg,1);
	obj.K1 = zeros(obj.nele,1);
    obj.Ex_num = zeros(obj.nedg,1);
    obj.Ey_num = zeros(obj.nedg,1);
    obj.Jx_num = zeros(obj.nedg,1);
    obj.Jy_num = zeros(obj.nedg,1);
    obj.H_num = zeros(obj.nele,1);
    obj.K_num = zeros(obj.nele,1);
end

