function [lb, ub] = CombineBounds(LB,UB)

num_state = size(LB,1);
if num_state ==1
    if any(isnan(LB))
        lb = nan;
        ub = nan;
    else
        lb = max(LB,[],'all');
        ub = min(UB,[],'all');
        if lb > ub
            lb = nan;
            ub = nan;
        end
    end
else
    lb = zeros(num_state,1);
    ub = lb;
    ind_nan_lb1 = find(isnan(LB(:,1)));
    ind_nan_lb2 = find(isnan(LB(:,2)));
    ind_nan_lb = union(ind_nan_lb1,ind_nan_lb2);
    ind_nan_ub1 = find(isnan(UB(:,1)));
    ind_nan_ub2 = find(isnan(UB(:,2)));
    ind_nan_ub = union(ind_nan_ub1,ind_nan_ub2);
    ind_nan = union(ind_nan_lb,ind_nan_ub);
    LB(ind_nan,:) = nan;
    UB(ind_nan,:) = nan;
    
    lb = max(LB,[],2);
    ub = min(UB,[],2);
    ind_empty = find(lb>ub);
    lb(ind_empty) = nan;
    ub(ind_empty) = nan;

    
    
end

end

