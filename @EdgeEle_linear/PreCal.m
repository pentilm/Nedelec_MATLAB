function PreCal( obj )
obj.Mref = [2 0 -1 0;0 2 0 -1;-1 0 2 0;0 -1 0 2];   % scale: *dx*dy/6
obj.Mcurl = [obj.dx;obj.dy;obj.dx;obj.dy];  %  should be changed for non-uniform and triangle element
obj.ele_area(:,1) = obj.dx*obj.dy;   % for triangle or non-uniform mesh, should loop for every element
end

