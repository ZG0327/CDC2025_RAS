function [us] = proj_feasible_ctrl(us,a1,b1,a2,b2,umax)
% us: 2-by-1, a:1-by-2
dim = size(us); 
% If any bs are NaN, this means already in infeasible region, therefore ACS
% should be NAN
if isnan(b1) || isnan(b2) 
    us = NaN(dim);

% If both a1 a2 = 0, there's the ACS is just U
elseif all(a1 == 0) && all(a2 == 0)
    us = us; 

% If a1 = 0 and a2 ~= 0, project to a2  
elseif all(a1 == 0) && any(a2 ~= 0)
    us = proj_ACS(us,a2,b2,umax);

% If a2 = 0 and a1 ~= 0, project to a1  
elseif all(a2 == 0) && any(a1 ~= 0)
    us = proj_ACS(us,a1,b1,umax);

% If both a1 a2 ~= 0, directly solve for [a1';a2']u=[b1;b2], and check if
% stays in U
elseif any(a1 ~= 0) && any(a2 ~= 0)
    A = [a1;a2];
    B = [b1;b2];
    sol = A\B;
    if any(abs(sol)-1e-4>umax)
        us = [NaN;NaN];
    else 
        us = sol;
    end

end

end

function u = proj_ACS(u,a,b,umax)
    % first see if U intersects with a,b
    vertices = [umax,umax; umax,-umax; -umax,umax; -umax,-umax];
    if all(vertices*a'>b)
        [~,ind] = min(vertices*a'-b);
        u = vertices(ind,:)';
    elseif all(vertices*a'-1e-4<=b)
        u = u;
    else
        % a'ubar > b, a'u = b , u-ubar = alpha a --> b - a'ubar = alpha a'a
        alpha = (b - a*u)/(a*a');
        u_temp = u+ alpha*a' ;

        if any(abs(u_temp)>umax)
            ind = find(abs(u_temp)>umax);
            dist_ind = abs(u_temp(ind)) - umax;
            P = [0,-1;1,0];
            abar = a*P;
            alphabar = dist_ind/abar(ind);
            u = u_temp -sign(u_temp(ind))*abar'*alphabar;
        else
            u = u_temp;
        end
    end

end