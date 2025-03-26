function dVdt = find_time_deriv(g,V,x,dt,t_reach)

% Remember to flip the data!!!!!!

if nargin<5
    error('specify the time vector')
end

clns = repmat({':'}, 1, g.dim);
t = size(V,g.dim+1);

if isvector(x)
    dVdt = -(eval_u(g,V(clns{:},t_reach),x)-eval_u(g,V(clns{:},t_reach+1),x))/dt;
else
    dVdt_all = zeros(prod(g.N, 'all'),1);
    dVdt = zeros(prod(g.N, 'all'),1);
    for i = 1 : t-1
        dVdt_all = -(eval_u(g,V(clns{:},i),x)-eval_u(g,V(clns{:},i+1),x))/dt;
        ind = find(t_reach == i);
        dVdt(ind) = dVdt_all(ind);
    end
    % dVdt = dVdt_all(:,t_reach);

end


end
