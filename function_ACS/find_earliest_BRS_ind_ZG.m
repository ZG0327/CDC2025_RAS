function tEarliest = find_earliest_BRS_ind_ZG(g, data, x, upper, lower)
% tEarliest = find_earliest_BRS_ind(g, data, x, upper, lower)
%     Determine the earliest time that the current state is in the reachable set
% 
% Inputs:
%     g, data - grid and value function representing reachable set
%     x       - state of interest
%     upper, lower - upper and lower indices of the search range
%
% Output:
%     tEarliest - earliest time index that x is in the reachable set
num_state = size(x,1);
small = 1e-4;
clns = repmat({':'}, 1, g.dim);

if num_state == 1
    if nargin < 4
        upper = size(data, g.dim+1);
        lower = 1;
    end
    % Binary search
   

    while upper > lower
        tEarliest = ceil((upper + lower)/2);
        valueAtX = eval_u(g, data(clns{:}, tEarliest), x);

        if valueAtX < small
            % point is in reachable set; eliminate all lower indices
            lower = tEarliest;
        else
            % too late
            upper = tEarliest - 1;
        end
    end

tEarliest = size(data, g.dim+1) - upper + 1;
% tEarliest = upper ;

else
    upper = repmat(size(data, g.dim+1),num_state,1);
    lower = ones(num_state,1);
    tEarliest = upper;

    ind_finished = [];
    ind_all = [1:num_state]';
    t = 1;

    while t <= size(data, g.dim+1) || length(ind_finished) == num_state
        % ind_remaining = setdiff(ind_all,ind_finished);

        % tEarliest = ceil((upper +lower)/2);

        valueAtX = eval_u(g, data(clns{:}, t), x);
        
        ind_reached = find(valueAtX<small);
        ind_reached = setdiff(ind_reached,ind_finished);
        tEarliest(ind_reached) = t;
        ind_finished = union(ind_finished,ind_reached);

        t = t+1;
    end    
    % tEarliest = size(data, g.dim+1) - tEarliest +1;
end


end