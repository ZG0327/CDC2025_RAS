function [lb,ub] = findUpLowBound(a,b,l_bound,u_bound)
% a & b defines the hyperplane, we use them to find the lb and ub

if ~isequal(class(a),class(b))
    error('a and b must have the same date type')
end

small = 0;

if isnan(a) 
    lb = nan;
    ub = nan;


elseif ~isscalar(a)
    u = b./a;
    ind_neg = find(a < -small);
    ind_pos = find(a > small);
    ind_0 = find( abs(a) <= small);
    ind_N = find(isnan(a));

    ub = nan(size(a));
    lb = ub;

    ub(ind_0) = u_bound;
    lb(ind_0) = l_bound;

    ub(ind_pos) = max(l_bound,min(u_bound,u(ind_pos)));
    lb(ind_pos) = l_bound;

    ub(ind_neg) = u_bound;
    lb(ind_neg) = min(u_bound,max(l_bound,u(ind_neg)));

    ub(ind_N) = nan;
    lb(ind_N) = nan;
elseif isscalar(a)
    u = b/a;
    if abs(a) <= small
        lb = l_bound;
        ub = u_bound;
    elseif a < -small 
        ub = u_bound;
        lb = min(u_bound,max(l_bound,u));
    else 
        lb = l_bound;
        ub = max(l_bound,min(u_bound,u));
    end
end


end