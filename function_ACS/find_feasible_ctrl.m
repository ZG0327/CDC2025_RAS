function [us] = find_feasible_ctrl(us,a1,b1,a2,b2,umax,umin)

if ~isscalar(b1)

% Get rid of those not in reach/avoid BRS.
    ind_reach_nan = find(isnan(b1));
    ind_avoid_nan = find(isnan(b2));
    ind_nan = union(ind_reach_nan,ind_avoid_nan);
    us(ind_nan,:) = nan;
% Note a1 b1 are for reach, a2 b2 are for avoid
res = sum(a1.*us,2) - b1;
ind_infeas = find(res > 0);
ind_infeas = setdiff(ind_infeas,ind_nan);

% First, make all control feasible for reach problem
alpha = -res(ind_infeas)./(a1(ind_infeas,1).*a1(ind_infeas,1)+...
    a1(ind_infeas,2).*a1(ind_infeas,2));

us(ind_infeas,:) = us(ind_infeas,:)+alpha.*a1(ind_infeas,:);
P = [0,-1;1,0];
abar = a1*P;
ind_out1 = find(abs(us(:,1))>1);
dist_ind1 = abs(us(ind_out1,1)) - 1;
alphabar = dist_ind1./abar(ind_out1,1);
us(ind_out1,:) = us(ind_out1,:) - ...
    sign(us(ind_out1,1)).*abar(ind_out1,:).*alphabar;

ind_out2 = find(abs(us(:,2))>1);
dist_ind2 = abs(us(ind_out2,2)) - 1;
alphabar = dist_ind2./abar(ind_out2,2);
us(ind_out2,:) = us(ind_out2,:) - ...
    sign(us(ind_out2,2)).*abar(ind_out2,:).*alphabar;

% Now, make all control feasible both on reach and avoid (if possible)
res_a = sum(a2.*us,2) - b2;
ind_infeas_a = find(res_a > 0);
ind_infeas_a = setdiff(ind_infeas_a,ind_nan);
% first check if reach constraints exists, than simply project to avoid
% constarints
ind_r_active = find(b1~=0 & ~isnan(b1));
ind_temp = setdiff(ind_infeas_a,ind_r_active);
alpha2 = -res(ind_temp)./(a2(ind_temp,1).*a2(ind_temp,1)+...
    a2(ind_temp,2).*a2(ind_temp,2));

us(ind_temp,:) = us(ind_temp,:)+alpha2.*a2(ind_temp,:);
abar2 = a2*P;
ind_out3 = find(abs(us(:,1))>1);
dist_ind3 = abs(us(ind_out3,1)) - 1;
alphabar2 = dist_ind3./abar2(ind_out3,1);
us(ind_out3,:) = us(ind_out3,:) - ...
    sign(us(ind_out3,1)).*abar2(ind_out3,:).*alphabar2;

ind_out4 = find(abs(us(:,2))>1);
dist_ind4 = abs(us(ind_out4,2)) - 1;
alphabar2 = dist_ind4./abar2(ind_out4,2);
us(ind_out4,:) = us(ind_out4,:) - ...
    sign(us(ind_out4,2)).*abar2(ind_out4,:).*alphabar2;

% if both exists, solve [a1';a2']u = [b1;b2] 
ind_both = setdiff(ind_infeas_a,ind_temp);
if ~isempty(ind_both)
    for i = ind_both'
        A = [a1(i,:); a2(i,:)];
        b = [b1(i);b2(i)];
        sol = A\b;
        if any(isnan(sol)) | any(abs(sol)>1)
            us(i,:) = [NaN,NaN];
        else
            us(i,:) = sol';
        end
    end
end

else 
    % res = sum(a1*us,2) - b1;
    if b1 == 0 && b2 == 0
       us = us;
    elseif b1 ~= 0 && b2 == 0
        alpha = (b1 - a1*us)/(a1*a1');
        u_online = us+ alpha*a1';
        P = [0,-1;1,0];
        abar = a1*P;
        if any(abs(u_online)>1)
        ind = find(abs(u_online)>1);
        dist_ind = abs(u_online(ind)) - 1;
        alphabar = dist_ind/abar(ind);
        us = u_online -sign(u_online(ind))*abar'*alphabar
        else
            us = u_online;
        end
    elseif b1 == 0 && b2 ~= 0
        alpha = (b2 - a2*us)/(a2*a2');
        u_online = us+ alpha*a2';
        P = [0,-1;1,0];
        abar = a2*P;
        ind = find(abs(u_online)>1);
        dist_ind = abs(u_online(ind)) - 1;
        alphabar = dist_ind/abar(ind);
        us = u_online -sign(u_online(ind))*abar'*alphabar;

    elseif b1 ~= 0 && b2 ~= 0
        A = [a1; a2];
        b = [b1;b2];
        sol = A\b;
        if any(isnan(sol)) | any(abs(sol)>1)
            us = [NaN;NaN];
        else
            us = sol;
        end
    end
    %     % res2 = sum(a2*u_fea,2) - b2;
    %     alpha2 = (b2 - a2*us)/(a2*a2');
    %     u_online2 = us+ alpha2*a2';
    %     P = [0,-1;1,0];
    %     abar2 = a2*P;
    %     ind2 = find(abs(u_online2)>1);
    %     dist_ind2 = abs(u_online2(ind2)) - 1;
    %     alphabar2 = dist_ind2/abar2(ind2);
    %     us = u_online2 -sign(u_online2(ind2))*abar2*alphabar2;
    % 
    % 
    % else
    %     A = [a1; a2];
    %     b = [b1;b2];
    %     sol = A\b;
    %     if any(isnan(sol)) | any(abs(sol)>1)
    %         us = [NaN,NaN];
    %     else
    %         us = sol';
    %     end
    % end

end