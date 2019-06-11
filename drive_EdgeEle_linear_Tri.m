% clc;
% clear;
maxwell = EdgeEle_linear_Tri([5,5,10]);
maxwell.Solve();
disp(maxwell.err_l2);
%disp(maxwell.err_inf);
