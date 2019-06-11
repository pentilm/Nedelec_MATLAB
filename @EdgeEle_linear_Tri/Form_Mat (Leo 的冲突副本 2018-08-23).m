function Form_Mat( obj )
% form matM and matBM, just run one time
% matM;   % global mass matrix
% matBM;  % global (1,curl phi)
% RHS;    % global right hand side vector, time independent

for i = 1:obj.nele     % loop for all elements
%             (x1,y1)
%               / \
%              /   \
%             /     \
%        (x2,y2)---(x3,y3)    
%
   x1 = obj.no2xy(1,obj.el2no(1,i));  y1 = obj.no2xy(2,obj.el2no(1,i));
   x2 = obj.no2xy(1,obj.el2no(2,i));  y2 = obj.no2xy(2,obj.el2no(2,i));
   x3 = obj.no2xy(1,obj.el2no(3,i));  y3 = obj.no2xy(2,obj.el2no(3,i));
% the gradient vectors scaled by 2|K|, see the document
v = zeros(2,3);    % v = [v1,v2,v3] where each v_i is a column vector
v(:,1) = [y2-y3;x3-x2];
v(:,2) = [y3-y1;x1-x3];
v(:,3) = [y1-y2;x2-x1];
for j = 1:obj.side_ele     % loop for all edges in the element
ed1 = obj.el2ed(j,i);
%% assemble for mass and stiffness matrix
% stiffness matrix
j1 = mod(j,obj.side_ele)+1;
obj.matBM(ed1,i) = obj.matBM(ed1,i) ... 
    + obj.edori(i,j)*obj.Mcurl(j)/obj.ele_area(i)*det([v(:,j),v(:,j1)]);    % be careful for the Mcurl
% mass matrix
    % for k == j, diagonal entries. notice edori are repeated for the same
    % edge, so their product must be 1. Ignored in the code.
    obj.matM(ed1,ed1) = obj.matM(ed1,ed1) ...
        + obj.Mref(j,j)/obj.ele_area(i)*(dot(v(:,j),v(:,j))+dot(v(:,j1),v(:,j1))-dot(v(:,j),v(:,j1)));
    
    % Since symmetry of the matrix has been used, just consider k > j
    for k = (j+1):obj.side_ele    % symmetry of mass matrix has been used here
        k1 = mod(k,obj.side_ele)+1;
        % to write down the formula, should reform the point indeces.
        % same notations to the document
        ej = intersect([j,j1],[k,k1]); % the common point
        ei = setdiff([j,j1],ej);   % the other point for the 1st edge
        ek = setdiff([k,k1],ej);   % the other point for the 2nd edge
        ed2 = obj.el2ed(k,i);
           
        % change here for non-uniform or triangle element
        % be careful for the Mcurl
        obj.matM(ed1,ed2) = obj.matM(ed1,ed2) + obj.edori(i,j)*obj.edori(i,k)*obj.Mref(j,k)/obj.ele_area(i) ... 
            *( dot(v(:,ei),v(:,ej))+dot(v(:,ej),v(:,ek))-dot(v(:,ej),v(:,ej))-2*dot(v(:,ei),v(:,ek)) );
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

