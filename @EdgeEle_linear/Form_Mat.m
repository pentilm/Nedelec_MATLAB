function Form_Mat( obj )
% form matM and matBM, just run one time
% matM;   % global mass matrix
% matBM;  % global (curl phi,1)
% RHS;    % global right hand side vector, time independent

for i = 1:obj.nele     % loop for all elements
% coordinates of this element from 1st node & 3rd node
%
%    (xae,ybe)----(xbe,ybe)
%       |           |
%       |           |
%    (xae,yae)----(xbe,yae)
%
   xae = obj.no2xy(1,obj.el2no(1,i));  xbe = obj.no2xy(1,obj.el2no(3,i));
   yae = obj.no2xy(2,obj.el2no(1,i));  ybe = obj.no2xy(2,obj.el2no(3,i));
obj.ele_area(i) = (xbe - xae)*(ybe - yae);  % compute the area of the element
for j = 1:obj.side_ele     % loop for all edges in the element
    ed1 = obj.el2ed(j,i);
%% assemble for mass and stiffness matrix
% stiffness matrix
obj.matBM(ed1,i) = obj.matBM(ed1,i) + obj.edori(i,j)*obj.Mcurl(j);    % be careful for the Mcurl
% mass matrix
    for k=j:obj.side_ele    % symmetry of mass matrix has been used here!
        ed2 = obj.el2ed(k,i);
        
        % change here for non-uniform or triangle element
        % be careful for the Mcurl
        obj.matM(ed1,ed2) = obj.matM(ed1,ed2) + obj.edori(i,j)*obj.edori(i,k)*obj.Mref(j,k)*obj.dx*obj.dy/6;   
        obj.matM(ed2,ed1) = obj.matM(ed1,ed2);     % symmetry!
    end

end
end
obj.matC = spdiags(obj.ele_area,0,obj.nele,obj.nele);   % the coefficient matrix for H
obj.matC_inv = spdiags(1./obj.ele_area,0,obj.nele,obj.nele);   % inverse of matC

% for zero Dirichlet B.C. take out interior elements
obj.matM = obj.matM(obj.eint,obj.eint);
obj.matBM = obj.matBM(obj.eint,:);

end

