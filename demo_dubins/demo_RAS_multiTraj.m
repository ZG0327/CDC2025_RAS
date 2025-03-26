close all; clear; clc

%% Import data
g_grid = importdata('grid.mat');
clns = repmat({':'}, 1, g_grid.dim);

u_lb = -pi;
u_ub = pi;


V_avoid_full = importdata('V_avoid.mat');
V_clvf = importdata('V_clvf.mat');
V_reach = importdata('V_reach.mat');
V_R = flip(V_reach,g_grid.dim+1);
V_reach2 = importdata('V_reach2.mat');
V_R2 = flip(V_reach2,g_grid.dim+1);
V_reach_avoid = importdata('V_ori_RA.mat');

V_avoid = V_avoid_full(:,:,:,end);

obs = V_avoid_full(clns{:},1);
goal = V_reach(clns{:},1);
goal2 = V_reach2(clns{:},1);

f_dyn = cell(3,1);
f_dyn{1} = cos(ones(g_grid.N'));
f_dyn{2} = sin(ones(g_grid.N'));
f_dyn{3} = zeros(g_grid.N');

g_dyn = cell(3,1);
g_dyn{1} = zeros(g_grid.N');
g_dyn{2} = zeros(g_grid.N');
g_dyn{3} = ones(g_grid.N');

gd_dyn = cell(3,2);
gd_dyn{1,1} = ones(g_grid.N');
gd_dyn{2,1} = zeros(g_grid.N');
gd_dyn{3,1} = zeros(g_grid.N');

gd_dyn{1,2} = zeros(g_grid.N');
gd_dyn{2,2} = ones(g_grid.N');
gd_dyn{3,2} = zeros(g_grid.N');
dRange = [-0.2,0.2];

%% Online simulation
dt = 0.1;
sim_t = [0:dt:18];

x0 = [-2.7;-0.4;1.9];

check_V = eval_u(g_grid,V_reach(clns{:},51),x0);
x = nan(3,length(sim_t));
u = nan(1,length(sim_t));
d = nan(2,length(sim_t));
x(:,1) = x0;
t = 0;

%%

% figure
% hold off
tic
for i = 1 : length(sim_t)
    i
    if i <=50
        t_reach = find_earliest_BRS_ind(g_grid, V_R, x(:,i), 51, 1);
        if t_reach>i
            x(:,i+1) = x(:,i);
        else
            dVdt_R = -(eval_u(g_grid,V_R(clns{:},t_reach),x(:,i))-...
                eval_u(g_grid,V_R(clns{:},t_reach+1),x(:,i)))/dt;

            exrReach.dVdt = dVdt_R;
            exrReach.x = x(:,i);
            exrReach.threshold = [-0.1,0.01];
            exrReach.dRange = dRange;
            [a_r,b_r] = findACS_D(V_R(clns{:},i),g_grid,f_dyn,g_dyn,gd_dyn,exrReach);
            [lb_r,ub_r] = findUpLowBound(a_r,b_r,u_lb,u_ub);
            % exrAvoid.dVdt = dVdt_A;
            exrAvoid.x = x(:,i);
            exrAvoid.dRange = dRange;
            [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrAvoid);
            [lb_a,ub_a] = findUpLowBound(a_a,b_a,u_lb,u_ub);
            [lb,ub] = CombineBounds([lb_r,lb_a],[ub_r,ub_a]);
            u(i) = (lb+ub)/2 + (rand(1)-0.5)*(ub-lb);
            check_V = eval_u(g_grid,V_R(clns{:},i),x(:,i))
            if isnan(u(i))
                break
            end
            d(:,i) = 2*(rand-0.5)/5;
            [~, x_temp] = ode45(@(t, s) Dubins_dyn(t, s, [d(:,i);u(i)]), [t t+dt], x(:,i));
            x(:,i+1) = x_temp(end,:);
            t = t+dt;
        end
    elseif i > 50 && i<= 130
        t_reach = find_earliest_BRS_ind(g_grid, V_R2, x(:,i), 81, 1)
        if t_reach>i-50
            x(:,i+1) = x(:,i);
        else
            dVdt_R = -(eval_u(g_grid,V_R2(clns{:},t_reach),x(:,i))-...
                eval_u(g_grid,V_R2(clns{:},t_reach+1),x(:,i)))/dt;

            exrReach.dVdt = dVdt_R;
            exrReach.x = x(:,i);
            exrReach.threshold = [-0.1,0.02];
            exrReach.dRange = dRange;
            [a_r,b_r] = findACS_D(V_R2(clns{:},i-50),g_grid,f_dyn,g_dyn,gd_dyn,exrReach);
            [lb_r,ub_r] = findUpLowBound(a_r,b_r,u_lb,u_ub);
            % exrAvoid.dVdt = dVdt_A;
            exrAvoid.x = x(:,i);
            exrAvoid.dRange = dRange;
            [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrAvoid);
            [lb_a,ub_a] = findUpLowBound(a_a,b_a,u_lb,u_ub);
            [lb,ub] = CombineBounds([lb_r,lb_a],[ub_r,ub_a]);
            u(i) = (lb+ub)/2 + (rand(1)-0.5)*(ub-lb);
            check_V = eval_u(g_grid,V_R2(clns{:},i-50),x(:,i))
            if isnan(u(i))
                break
            end
            d(:,i) = 2*(rand-0.5)/5;
            [~, x_temp] = ode45(@(t, s) Dubins_dyn(t, s, [d(:,i);u(i)]), [t t+dt], x(:,i));            
            x(:,i+1) = x_temp(end,:);
            t = t+dt;
        end
    else
        exrCLF.x = x(:,i);
        exrCLF.gamma = 0.1;
        exrA.gamma = 0;
        exrA.x = x(:,i);
        exrA.dRange = dRange;
        exrCLF.dRange = dRange;
        [a_clvf,b_clvf] = findACS_D(V_clvf,g_grid,f_dyn,g_dyn,gd_dyn,exrCLF);
        [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrA);

        [lb_clvf,ub_clvf] = findUpLowBound(a_clvf,b_clvf,u_lb,u_ub);
        [lb_a,ub_a] = findUpLowBound(a_a,b_a,u_lb,u_ub);

        [lb,ub] = CombineBounds([lb_clvf,lb_a],[ub_clvf,ub_a]);

        d(:,i) = 2*(rand-0.5)/5;
        u(i) = (lb+ub)/2 + (rand(1)-0.5)*(ub-lb);
        check_V = eval_u(g_grid,V_clvf,x(:,i))
        % if check_V <= -0.03
        %     break
        % end
        [~, x_temp] = ode45(@(t, s) Dubins_dyn(t, s, [d(:,i);u(i)]), [t t+dt], x(:,i));
        x(:,i+1) = x_temp(end,:);
        t = t+dt;

    end


end

toc


%%

figure
hold on
plot3(x(1,:),x(2,:),x(3,:),'k')
% visSetIm(g_grid,goal,'m',0);
visSetIm(g_grid,goal2,'m',0);

visSetIm(g_grid,obs,'r',0);
visSetIm(g_grid,V_avoid,'r',0);
% visSetIm(g_grid,V_clvf,'k',0.5)
% visSetIm(g_grid,V_reach2(clns{:},51),'b',0);
% visSetIm(g_grid,V_reach(clns{:},51),'b',0);

% visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0);

% plot(g_grid.xs{1}(ind_infeasible),g_grid.xs{2}(ind_infeasible),'k.')
% plot(g_grid.xs{1}(ind_nan),g_grid.xs{2}(ind_nan),'r.')


%%
function dydt = Dubins_dyn(t,y,u)

dydt = [cos(y(3))+u(1); sin(y(3))+u(2) ; u(3)];
end
