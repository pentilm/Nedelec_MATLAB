function Partition2D( obj )

% define the corner nodes
node = [obj.ax,obj.ay;
        obj.bx,obj.ay;
        obj.bx,obj.by;
        obj.ax,obj.by];

% test density function for non-uniform mesh
% hdata.fun = @(x,y)(0.01 + 0.1*sqrt( (x-0.25).^2+(y-0.75).^2 ))*10;

% generate triangle mesh
hdata.hmax = sqrt((obj.bx-obj.ax)*(obj.by-obj.ay)/(obj.Nx-1)/(obj.Ny-1));   % approximate the mesh size
[p,t] = mesh2d(node,[],hdata);
[e,te,e2t,bnd] = connectivity ( p, t );
obj.bnd = bnd;
% get the basic information of the mesh
obj.npt = size(p,1);
obj.nele = size(t,1);
obj.nedg = size(e,1);
% initialize all the arries
Allo(obj);
% get area for every element
obj.ele_area(:,1) = triarea(p,t);

% generate a triangle mesh
obj.no2xy = p';

% initialize the mid points array
% triangle partition, find the centroid of every triangle
c = centroid_mesh ( p, t );
obj.px_mid = c(:,1)';
obj.py_mid = c(:,2)';

% 3 nodes (counterclockwise) for each element!
obj.el2no = t';

% the total number of edges including boundary edges
obj.ed2no = e;
obj.el2ed = te';

% Indicators: 0 for interior edges; 1 for boundary edges.
% notice all the entries of ed_id has been initialized by 1
for i = 1:obj.nedg
    if (e2t(i,2)==0)
       obj.ed_id(i) = 1; 
    end
end

obj.eint = find(obj.ed_id == 0); % get labels for all interior edges
obj.ienum = length(obj.eint);    % total number of interior edges

% total unknown = # of H (elements) + # of E (inner edges)
obj.num_unk = obj.nele + obj.ienum;

% compare reference element edge directions vs global 
% edge directions to get the orientations for all edges
% criteria: In array 'e' (got by connectivity), we already got a global
% edge direction. Now compare each edge in every element to this global
% direction counter-clockwisely.

for i = 1:obj.nele
    for j = 1:obj.side_ele
        % given j-th edge of i-th triangle, find its edge number
        % to get 1->1, 2->2, 3->3, 4->1, use: mod(n-1,3)+1
        edg_num = (ones(obj.nedg,1)*sort(t(i,[j,mod(j,obj.side_ele)+1]),2) == e);
        % if 1st and 2nd elements of edg_num are all 1, then this is the
        % target edge in e
        if t(i,j) ~= e(edg_num(:,1) == 1 & edg_num(:,2) == 1 ,1)   % if the order in t is same to the order in e
            % if 1st element of t does not equal 2nd element in e, then change the
            % sign
            obj.edori(i,j) = -1;
        end
    end                 
end
end

