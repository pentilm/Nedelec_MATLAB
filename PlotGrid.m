function PlotGrid( x, y, labelflag )
% Plot grid which the data given by x and y
% x and y must be 1d vectors, no matter row or column
% labelflag = 1: title the x and y labels
if (nargin<3)
        labelflag = 1;
end
plot(x,meshgrid(y,x),'k');
hold on;
plot(meshgrid(x,y),y,'k');
axis([min(x),max(x),min(y),max(y)]);
if labelflag
   xlabel('x');
   ylabel('y');
end
end

