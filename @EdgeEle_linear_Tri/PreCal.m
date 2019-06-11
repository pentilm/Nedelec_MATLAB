function PreCal( obj )
% all the elements here must be scaled when use. see the document
obj.Mref = [1/24,1/48,1/48;1/48,1/24,1/48;1/48,1/48,1/24];
obj.Mcurl = 0.5*ones(3,1);
end

