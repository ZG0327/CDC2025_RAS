function [a,b] = findACS_D(V,grid,f,g,gd,exrArg)
% V: value function
% grid: grid
% f,g: dynamics of the control affine system. n-by-1 (n-by-m) cell, where n
%       is the number of state, m is the number of control
% exrArg.
%       timeVarying: ('boolean') 1: find U(x,t), 0: find U(x).
%       gamma: specifally for CLVF
%       dVdt
%       dRange = [dmin, dmax]

%
    if isfield(exrArg,'gamma')
        gamma = exrArg.gamma;
    else
        gamma = 0;
    end

    if isfield(exrArg,'dRange')
        dRange = exrArg.dRange;
    else
        error('specify dRange');
    end
    
    if isfield(exrArg,'x')
        x = exrArg.x;
    else
        x = [];
    end

    if isfield(exrArg,'dVdt')
        dVdt = exrArg.dVdt;
    else
        if isvector(x)  
            dVdt = 0;
        elseif ~isvector(x)
            dVdt = zeros(size(V));
        end
    end

    if isfield(exrArg,'threshold')
        small = exrArg.threshold;
    else
        small= [-0.01,0.01];
    end


    %%
    if isvector(x)  
        % finds a/b for one state with Time invariant functions
        [a,b]  = find_State_ACS(V,grid,f,g,gd,dRange,x,gamma,dVdt,small);
    elseif ~isvector(x) 
        % finds a/b for all grid points with Time invariant functions
        [a,b] = find_grid_ACS(V,grid,f,g,gd,dRange,x,gamma,dVdt,small);
    end
end

function [a,b] = find_grid_ACS(V,grid,f,g,gd,dRange,x,gamma,dVdt,small)

% if nargin <6
%     dVdt = 0;
% end

Deriv = computeGradients(grid,V);
dim = size(Deriv,1);
dVdt = reshape(dVdt,prod(grid.N, 'all'),[]);
dim_gd = size(gd);
% xs = zeros(prod(grid.N, 'all'),2);
% xs(:,1) = reshape(grid.xs{1}, 1, []);
% xs(:,2) = reshape(grid.xs{2}, 1, []);
xs = x;
% V_xs = reshape(V, 1, [])';
Vs = eval_u(grid,V,xs);

deriv_xs = cell(size(Deriv));
f_xs = cell(size(f));
g_xs = cell(size(g));
gd_xs = cell(size(gd));

LfV = zeros(prod(grid.N, 'all'),1);
LgV = LfV;
LgdV = zeros(prod(grid.N, 'all'),dim_gd(2)); 
d = zeros(prod(grid.N, 'all'),1);

for i = 1:dim
    deriv_xs{i} = eval_u(grid,Deriv{i},xs);
    f_xs{i} = eval_u(grid,f{i},xs);
    g_xs{i} = eval_u(grid,g{i},xs);
    for j = 1:dim_gd(2)
    gd_xs{i,j} = eval_u(grid,gd{i,j},xs);
    LgdV(:,j) = LgdV(:,j) + deriv_xs{i}.*gd_xs{i,j};
    end
    LfV = LfV + deriv_xs{i}.*f_xs{i};
    LgV = LgV + deriv_xs{i}.*g_xs{i};
    
end

% Compute 'a(x)' and 'b(x)'
a = nan(prod(grid.N, 'all'),1);
b = a;

% ind_d_pos = LgdV>0;
% ind_d_neg = LgdV<=0;
% 
% d(ind_d_pos) = dRange(2);
% d(ind_d_neg) = dRange(1);

d = (LgdV>0)*dRange(2)+(LgdV<=0)*dRange(1) ;
% lgdV_d = LgdV.*d;
lgdV_d = sum(LgdV.*d,2);

if gamma == 0
    ind_acs = find(Vs<small(2) & Vs>small(1));
    a(ind_acs) = LgV(ind_acs);
    b(ind_acs) = -dVdt(ind_acs)-LfV(ind_acs)-lgdV_d(ind_acs);

    ind_pos = find(Vs>=small(2));
    a(ind_pos) = nan;
    b(ind_pos) = nan;

    ind_neg = find(Vs<=small(1));
    a(ind_neg) = 0;
    b(ind_neg) = 0;
else
    a = LgV;
    b = -LfV-lgdV_d-gamma*Vs+0;
end

end


%%
function [a,b] = find_State_ACS(V,grid,f,g,gd,dRange,x,gamma,dVdt,small)           
    % 
    % if nargin<7
    %     dVdt = 0;
    % end
    % data process
    Deriv = computeGradients(grid,V);
    dim = size(Deriv,1);
    dim_gd = size(gd);
    % LfV = zeros(dim,1);
    % LgV = zeros(dim,1);
    % LgdV = zeros(dim,1);
    dVdx = zeros(dim,1);

    Vx = eval_u(grid,V,x);
    fx = eval_u(grid,f,x);
    gx = eval_u(grid,g,x);
    gdx = eval_u(grid,gd,x);
    for i = 1: dim_gd(2)
        gdx(:,i) = eval_u(grid,gd(:,i),x);
    end
    for i = 1:dim
        dVdx(i) = eval_u(grid,Deriv{i},x);
    end

    LfV_x = dVdx'*fx;
    LgV_x = dVdx'*gx;
    LgdV_x = dVdx'*gdx;
    % if LgdV_x > 0
    %     d = dRange(2);
    % else 
    %     d = dRange(1);
    % end
    d = (LgdV_x>=0)*dRange(2)+(LgdV_x<0)*dRange(1); 
    if gamma == 0
        if Vx<small(2) && Vx>small(1)
            a = LgV_x;
            b = -dVdt-LfV_x-LgdV_x*d';
        elseif Vx>=small(2)
            a = nan;
            b = nan;
        else
            a = 0;
            b = 0;
        end
    else
        a = LgV_x;
        b = -LfV_x-LgdV_x*d'-gamma*Vx+0.0001;
    end


end