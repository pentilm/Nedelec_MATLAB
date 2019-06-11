function Initialize( obj )
% intialize for rectangle region (uniform)
obj.dx = (obj.bx-obj.ax)/(obj.Nx-1);
obj.dy = (obj.by-obj.ay)/(obj.Ny-1);
obj.dt = (obj.bt-obj.at)/(obj.Nt-1);

end

