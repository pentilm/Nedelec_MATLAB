function Initialize( obj )
% notice for triangle mesh, there is no dx and dy
obj.dt = (obj.bt-obj.at)/(obj.Nt-1);
end

