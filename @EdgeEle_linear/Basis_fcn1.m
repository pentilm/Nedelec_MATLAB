function y0 = Basis_fcn1(x,y,nn,ne,obj)
% (x, # of element, # of edge)
% basis in H^1(curl) of rectangle element
% piecewise linear
% just return the non-zero component of basis function

% coordinates of this element from 1st node & 3rd node
xae = obj.no2xy(1,obj.el2no(1,nn));  xbe = obj.no2xy(1,obj.el2no(3,nn));
yae = obj.no2xy(2,obj.el2no(1,nn));  ybe = obj.no2xy(2,obj.el2no(3,nn));
switch ne
    case 1
        y0 = (ybe-y)/(ybe-yae);
    case 2
        y0 = (x-xae)/(xbe-xae);
    case 3
        y0 = -(y-yae)/(ybe-yae);
    otherwise
        y0 = -(xbe-x)/(xbe-xae);
end
end

