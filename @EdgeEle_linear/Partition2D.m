function Partition2D( obj )
% generate a rectangular mesh
for j=1:obj.Ny
    for i=1:obj.Nx
        ipt=obj.Nx*(j-1)+i;
        obj.no2xy(1,ipt)=obj.dx*(i-1);   
        obj.no2xy(2,ipt)=obj.dy*(j-1);
    end
end
% initialize the mid points array
% uniform rectangle partition
obj.px_mid = linspace(0.5*obj.dx,obj.bx-0.5*obj.dx,obj.Nx-1);
obj.py_mid = linspace(0.5*obj.dy,obj.by-0.5*obj.dy,obj.Ny-1);

% 4 nodes (counterclockwise) for each element!
obj.el2no=zeros(obj.side_ele,obj.nele);    
idx=1;
for i=1:(obj.Ny-1)	% number of columns to go through
    for j=1:(obj.Nx-1)
        obj.el2no(:,idx)=[j+(i-1)*obj.Nx; j+(i-1)*obj.Nx+1; j+obj.Nx*i+1; j+obj.Nx*i];
        idx = idx+1;
    end
end

% the total number of edges including boundary edges
for i=1:obj.nele
    for j=1:obj.side_ele		% for each element in each column
        if (j~=obj.side_ele)
            obj.edge_info((i-1)*obj.side_ele+j,:)=[obj.el2no(j,i), obj.el2no(j+1,i)];
        else
            obj.edge_info((i-1)*obj.side_ele+j,:)=[obj.el2no(j,i), obj.el2no(1,i)];
        end
    end
end
obj.edge_info=sort(obj.edge_info,2);                     
[obj.ed2no,~,obj.el2ed]=unique(obj.edge_info,'rows');   
obj.el2ed=reshape(obj.el2ed,obj.side_ele,obj.nele);  
% Indicators: 0 for interior edges; 1 for boundary edges.
for i=1:obj.nedg
    v1x = obj.no2xy(1,obj.ed2no(i,1));
    v1y = obj.no2xy(2,obj.ed2no(i,1));
    v2x = obj.no2xy(1,obj.ed2no(i,2));
    v2y = obj.no2xy(2,obj.ed2no(i,2));
    if (v1x==obj.ax && v2x==obj.ax) || (v1x==obj.bx && v2x==obj.bx) || (v1y==obj.ay && v2y==obj.ay) || (v1y==obj.by && v2y==obj.by)
        obj.ed_id(i)=1;
    end
end

obj.eint = find(obj.ed_id == 0); % get labels for all interior edges
obj.ienum = length(obj.eint);    % total number of interior edges

% total unknown = # of H (elements) + # of E (inner edges)
obj.num_unk = obj.nele + obj.ienum;

% compare reference element edge directions vs global 
% edge directions to get the orientations for all edges
% criteria: counter clockwise
% retangle:
%        3
%   ----<------
%   |          |
% 4 v          ^  2
%   |          |
%   ----->-----
%         1

for i = 1:obj.nele
 %   edori(i,:)=[1 1 -1 -1];
    for j = 1:obj.side_ele
%         edn = obj.el2ed(j,i);
%         n1=ed2no(edn,1); n2=ed2no(edn,2);
        if j < obj.side_ele
          m1 = obj.el2no(j,i); m2 = obj.el2no(j+1,i);
        else
          m1 = obj.el2no(j,i); m2 = obj.el2no(1,i);
        end
        if m1 > m2 
            obj.edori(i,j) = -1;
        end
    end                 
end



end

