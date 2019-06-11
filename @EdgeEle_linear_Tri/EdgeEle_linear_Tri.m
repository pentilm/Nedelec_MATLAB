classdef EdgeEle_linear_Tri < handle
    % Use linear edge element to solve Drude's model
    properties
        %%
        % number of sides of element
        % dim_ele = 4 means retangle
        % dim_ele = 3 means triangle
        side_ele;
        %% computational parameters
        Nx;Ny;Nt;   % number of partition points
        ax;bx;
        ay;by;
        at;bt;
%         dx;dy;
        dt;
        %% functions
        % exact solutions
        Ex_exa;Ey_exa;
        Jx_exa;Jy_exa;
        H_exa;K_exa;
        % RHS functions
        g1_1;g1_2;g2;g3_1;g3_2;g4;
        %% equation parameters
        e0;mu0;
        ge;we;
        gm;wm;
        %% mesh data
        npt;    % total number of points
        nele;   % total number of elements
        nedg;   % total number of edges
        
        % 2*n array where n is total nodal points number. global nodal
        % number -> coordinate. n-th column is n-th point's coordinate [x;y]
        no2xy;
        
        % 3*n array where n is total element number. (local node,
        % global element) -> global node. n-th column is n-th element
        el2no;
        
        % information of every edge: [start pt, end pt]
%         edge_info;
        
        % 3*n array where n is total element number. (local edge, global
        % element) -> global edge number
        el2ed;
        
        % n*2 array where n is total edge number. (global edge, 1 or
        % 2(number of endpoint)) -> global points number of the indicated
        % edge
        ed2no;
        
        % global edge number -> inner (0) or boundary edges (1)
        ed_id;
        
        % labels for all interior points
        eint;
        % total number of interior edges
        ienum;
        
        % get the direction w.r.t. reference element
        % (number of global element, number of local edge) -> direction
        edori;
        
        %  bnd : Nx1 logical array identifying boundary nodes. no2xy(:,i) is a boundary
        % node if BND(i)=TRUE.
        bnd;
        
        %% computational data
        Mref;   % unit local mass matrix, should scales when use (the matrix with mesh size 1)
        Mcurl;  % unit (1,curl N_i), should scales when use (the matrix with mesh size 1)
        matM;   % global mass matrix
        matBM;  % global (curl phi,1)
        matC;   % global (1,1)  (the coefficient matrix for H)
        matC_inv;   % inverse of matC
        ele_area; % area of every element
        num_unk;    % number of unkowns
        RHS;    % global right hand side vector, time independent
        RHS_tmp;    % global right hand side vector, time dependent. update on each time step.
        E0;H0;  % 1D vector of old coefficients of E and H
        E1;H1;  % 1D vector of new coefficients of E and H
        J0;K0;  % 1D vector of old coefficients of J and K
        J1;K1;  % 1D vector of old coefficients of J and K
        vec_g3;vec_g4;  % vectors used to update J and K
        err_inf;err_l2; % error of all components in L^\infty and L^2
        Ex_num;Ey_num;  % numerical solution values at final time
        Jx_num;Jy_num;
        H_num;K_num;
        px_mid;py_mid;  % mid points for every element. used for plot.
    end
    
    methods
        function obj = EdgeEle_linear_Tri(N)
            % triangle element
            obj.side_ele = 3;
            % parameters
            obj.Nx = N(1);
            obj.Ny = N(2);
            obj.Nt = N(3);
            obj.ax = 0; obj.bx = 1; % x direction
            obj.ay = 0; obj.by = 1; % y direction
            obj.at = 0; obj.bt = 1e-4; % t direction
            obj.e0 = 1; obj.mu0 = 1;
            obj.ge = 1; obj.we = 1;
            obj.gm = 1; obj.wm = 1;
            % exact solution setting
            obj.Ex_exa = @(x,y,t)sin(pi*y).*exp(-obj.ge*t);
            obj.Ey_exa = @(x,y,t)sin(pi*x).*exp(-obj.ge*t);
            obj.Jx_exa = @(x,y,t)sin(pi*y).*exp(-obj.ge*t)*obj.we^2.*t;
            obj.Jy_exa = @(x,y,t)sin(pi*x).*exp(-obj.ge*t)*obj.we^2.*t;
            obj.H_exa = @(x,y,t)1/pi*(cos(pi*x)-cos(pi*y)).*exp(-obj.ge*t).*(obj.we^2*t-obj.ge);
            obj.K_exa = @(x,y,t)1/pi*(cos(pi*x)-cos(pi*y)).*exp(-obj.ge*t).*obj.we^2*(0.5*obj.we^2*t.^2-obj.ge*t);
            obj.g1_1 = @(x,y,t)0;
            obj.g1_2 = @(x,y,t)0;
            obj.g2 = @(x,y,t)1/pi*(cos(pi*x)-cos(pi*y)).*exp(-obj.ge*t).*(-2*obj.ge*obj.we^2*t+obj.ge^2+obj.we^2+pi^2+0.5*obj.we^4*t.^2);
            obj.g3_1 = @(x,y,t)0;
            obj.g3_2 = @(x,y,t)0;
            obj.g4 = @(x,y,t)0;

%             obj.Ex_exa = @(x,y,t)x.*(1-x).*y.*(1-y);
%             obj.Ey_exa = @(x,y,t)0;
%             obj.Jx_exa = @(x,y,t)0;
%             obj.Jy_exa = @(x,y,t)0;
%             obj.H_exa = @(x,y,t)0;
%             obj.K_exa = @(x,y,t)0;
%             obj.g1_1 = @(x,y,t)0;
%             obj.g1_2 = @(x,y,t)0;
%             obj.g2 = @(x,y,t)-x.*(1-x).*(1-2*y);
%             obj.g3_1 = @(x,y,t)-x.*(1-x).*y.*(1-y).*obj.we^2;
%             obj.g3_2 = @(x,y,t)0;
%             obj.g4 = @(x,y,t)0;
            
%             obj.Ex_exa = @(x,y,t)0;
%             obj.Ey_exa = @(x,y,t)0;
%             obj.Jx_exa = @(x,y,t)0;
%             obj.Jy_exa = @(x,y,t)0;
%             obj.H_exa = @(x,y,t)x+y;
%             obj.K_exa = @(x,y,t)x-y;
%             obj.g1_1 = @(x,y,t)-1;
%             obj.g1_2 = @(x,y,t)1;
%             obj.g2 = @(x,y,t)x-y;
%             obj.g3_1 = @(x,y,t)0;
%             obj.g3_2 = @(x,y,t)0;
%             obj.g4 = @(x,y,t)obj.gm*(x-y)-obj.wm^2*(x+y);

            obj.err_inf = zeros(6,1);
            obj.err_l2 = zeros(6,1);
        end
        %% member methods
        % since for triangle mesh, cannot predict element number etc, so
        % after partition, allocate arraies
        Allo(obj);  
        Initialize(obj);
        Partition2D(obj);
        PreCal(obj);
        Form_Mat(obj);  % form mass matrix and RHS vector
        Init_Val(obj);
        Form_RHS(obj);  % form time independent RHS vector
        Form_RHS_tmp(obj,ct)    % form RHS contribution at t = ct. Just form the time dependent part.
%         Form_Bound(obj);    % form boundary contribution
        Marching(obj);
        UpdateJK(obj);
        Solve(obj);
        Post_Process(obj);
        Plot_Sol(obj);
        %% tool methods
        y = Basis_fcn1(x,y,nn,ne,vn,obj);    % (x, # of element, # of edge)
    end
    
end

