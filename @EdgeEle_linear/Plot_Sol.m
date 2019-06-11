function Plot_Sol( obj )
% Plot all solutions
[X,Y] =meshgrid(obj.px_mid,obj.py_mid);
figure(1);
mesh(X,Y,obj.Ex_num);
figure(2);
mesh(X,Y,obj.Ey_num);
figure(3);
mesh(X,Y,obj.Jx_num);
figure(4);
mesh(X,Y,obj.Jy_num);
figure(5);
mesh(X,Y,obj.H_num);
figure(6);
mesh(X,Y,obj.K_num);

end

