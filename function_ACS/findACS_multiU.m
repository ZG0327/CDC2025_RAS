function [a,b] = findACS_multiU(V,grid,f,g,exrArg)
% V: value function
% grid: grid
% f,g: dynamics of the control affine system. n-by-1 (n-by-m) cell, where n
%       is the number of state, m is the number of control
% exrArg.
%       timeVarying: ('boolean') 1: find U(x,t), 0: find U(x)
%
    if isfield(exrArg,'gamma')
        gamma = exrArg.gamma;
    else
        gamma = 0;
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

small= 0.01;

    %%
    if isvector(x)  
        % finds a/b for one state with Time invariant functions
        [a,b]  = find_State_ACS(V,grid,f,g,x,gamma,dVdt,small);
    elseif ~isvector(x) 
        % finds a/b for all grid points with Time invariant functions
        [a,b] = find_grid_ACS(V,grid,f,g,x,gamma,dVdt,small);
    end
end

function [a,b] = find_grid_ACS(V,grid,f,g,x,gamma,dVdt,small)

% if nargin <6
%     dVdt = 0;
% end

Deriv = computeGradients(grid,V);
dim = size(Deriv,1);
udim = size(g,2);
dVdt = reshape(dVdt,prod(grid.N, 'all'),[]);

% xs = zeros(prod(grid.N, 'all'),2);
% xs(:,1) = reshape(grid.xs{1}, 1, []);
% xs(:,2) = reshape(grid.xs{2}, 1, []);
xs = x;
% V_xs = reshape(V, 1, [])';
Vs = eval_u(grid,V,xs);

deriv_xs = cell(size(Deriv));
f_xs = cell(size(f));
g_xs = cell(size(g));

LfV = zeros(prod(grid.N, 'all'),1);
LgV = zeros(prod(grid.N, 'all'),udim);


for i = 1:dim
    deriv_xs{i} = eval_u(grid,Deriv{i},xs);
    f_xs{i} = eval_u(grid,f{i},xs);
    LfV = LfV + deriv_xs{i}.*f_xs{i};
    for j = 1:udim
        g_xs{i,j} = eval_u(grid,g{i,j},xs);
        LgV(:,i) = LgV(:,i) + deriv_xs{i}.*g_xs{i,j};
    end
end

% Compute 'a(x)' and 'b(x)'
a = nan(prod(grid.N, 'all'),udim);
b = nan(prod(grid.N, 'all'),1);

if gamma == 0
    ind_acs = find(Vs<small & Vs>-small);
    a(ind_acs,:) = LgV(ind_acs,:);
    b(ind_acs) = -dVdt(ind_acs)-LfV(ind_acs);

    ind_pos = find(Vs>=small);
    a(ind_pos,:) = nan;
    b(ind_pos) = nan;

    ind_neg = find(Vs<=-small);
    a(ind_neg,:) = 0;
    b(ind_neg) = 0;
else
    a = LgV;
    b = -LfV-gamma*Vs+0.001;
end

end


%%
function [a,b] = find_State_ACS(V,grid,f,g,x,gamma,dVdt,small)           
    % 
    % if nargin<7
    %     dVdt = 0;
    % end
    % data process
    Deriv = computeGradients(grid,V);
    dim = size(Deriv,1);
    udim = size(g,2);
    
    dVdx = zeros(dim,1);

    Vx = eval_u(grid,V,x);
    fx = eval_u(grid,f,x);
    gx = eval_u(grid,g,x);
    for i = 1:dim
        dVdx(i) = eval_u(grid,Deriv{i},x);
        for j = 1:udim
break
        end
    end

    LfV_x = dVdx'*fx;
    LgV_x = dVdx'*gx;

    if gamma == 0
        if Vx<small && Vx>-small
            a = LgV_x;
            b = -dVdt-LfV_x;
        elseif Vx>=small
            a = nan;
            b = nan;
        else
            a = 0;
            b = 0;
        end
    else
        a = LgV_x;
        b = -LfV_x-gamma*Vx+0.001;
    end


end