function Plot_Sol( obj )
% Plot all solutions
Ex_plot = obj.Ex_num;
Ey_plot = obj.Ey_num;
Jx_plot = obj.Jx_num;
Jy_plot = obj.Jy_num;
H_plot = zeros(obj.npt,1);
K_plot = zeros(obj.npt,1);
% auxiliary arrays used to count how many times to visit XX_vec vector
H_count = zeros(obj.npt,1);
K_count = zeros(obj.npt,1);
for i = 1:obj.nele
    for j = 1:obj.side_ele   % loop for all points of each element
        v = obj.el2no(j,i);
       % H
       H_plot(v) = H_plot(v) + obj.H_num(i);
       H_plot(v) = H_plot(v) + obj.H_num(i);
       H_count(v) = H_count(v) + 1;
       % K
       K_plot(v) = K_plot(v) + obj.K_num(i);
       K_plot(v) = K_plot(v) + obj.K_num(i);
       K_count(v) = K_count(v) + 1;
    end
end
% compute average and replace NaN to 0
H_plot = H_plot./H_count;
H_plot(isnan(H_plot)) = 0;
K_plot = K_plot./K_count;
K_plot(isnan(K_plot)) = 0;

% mesh the surface
figure;
trimesh(obj.el2no', obj.no2xy(1,:)' , obj.no2xy(2,:)' , Ex_plot);
figure;
trimesh(obj.el2no', obj.no2xy(1,:)' , obj.no2xy(2,:)' , Ey_plot);
figure;
trimesh(obj.el2no', obj.no2xy(1,:)' , obj.no2xy(2,:)' , Jx_plot);
figure;
trimesh(obj.el2no', obj.no2xy(1,:)' , obj.no2xy(2,:)' , Jy_plot);
figure;
trimesh(obj.el2no', obj.no2xy(1,:)' , obj.no2xy(2,:)' , H_plot);
figure;
trimesh(obj.el2no', obj.no2xy(1,:)' , obj.no2xy(2,:)' , K_plot);

end

