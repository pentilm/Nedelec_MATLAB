clc;
clear;
maxwell = EdgeEle_linear([5,5,10]*2);
maxwell.Solve();
disp(maxwell.err_l2);
disp(maxwell.err_inf);
