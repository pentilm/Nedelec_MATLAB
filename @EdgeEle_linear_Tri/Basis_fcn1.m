function y0 = Basis_fcn1(x,y,nn,ne,vn,obj)
% (x, # of element, # of edge, # of element of result)
% basis in H^1(curl) of triangle element
% piecewise linear

%             (x1,y1)
%               / \
%              /   \
%             /     \
%        (x2,y2)---(x3,y3)    
%
   x1 = obj.no2xy(1,obj.el2no(1,nn));  y1 = obj.no2xy(2,obj.el2no(1,nn));
   x2 = obj.no2xy(1,obj.el2no(2,nn));  y2 = obj.no2xy(2,obj.el2no(2,nn));
   x3 = obj.no2xy(1,obj.el2no(3,nn));  y3 = obj.no2xy(2,obj.el2no(3,nn));
switch ne
    case 1
        switch vn
            case 1
                y0 = (y3 - y)/obj.ele_area(nn)*0.5;
            case 2
                y0 = (x - x3)/obj.ele_area(nn)*0.5;
        end
    case 2
        switch vn
            case 1
                y0 = (y1 - y)/obj.ele_area(nn)*0.5;
            case 2
                y0 = (x - x1)/obj.ele_area(nn)*0.5;
        end
    otherwise
        switch vn
            case 1
                y0 = (y2 - y)/obj.ele_area(nn)*0.5;
            case 2
                y0 = (x - x2)/obj.ele_area(nn)*0.5;
        end
end
end

