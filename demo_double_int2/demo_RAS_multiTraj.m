close all; clear; clc

%% Import data
V_avoid_full = importdata('V_avoid.mat');
V_clvf_all = importdata('V_clvf.mat');
V_clvf = V_clvf_all(:,:,end);
V_reach = importdata('V_reach.mat');
V_R = flip(V_reach,3);
V_A = flip(V_avoid_full(:,:,1:end),3);
V_reach_avoid = importdata('V_ori_RA.mat');
V_reach2 = importdata('V_reach2.mat');
V_R2 = flip(V_reach2,3);
V_avoid = V_avoid_full(:,:,end);

g_grid = importdata('grid.mat');
obs = V_avoid_full(:,:,1);
goal = V_reach(:,:,1);
goal2 = V_reach2(:,:,1);

f_dyn = cell(2,1);
f_dyn{1} = g_grid.xs{2};
f_dyn{2} = zeros(g_grid.N');

g_dyn = cell(2,1);
g_dyn{1} = zeros(g_grid.N');
g_dyn{2} = ones(g_grid.N');

gd_dyn = cell(2,1);
gd_dyn{1} = ones(g_grid.N');
gd_dyn{2} = zeros(g_grid.N');

dRange = [-0.2,0.2];

%% Online simulation
dt = 0.01;
sim_t = [0:dt:2];

x0 = [-0.96;0.86];
x = nan(2,length(sim_t));
u = nan(1,length(sim_t));
d = nan(1,length(sim_t));
x(:,1) = x0;

t = 0;

%%

% figure
% hold off
tic
for i = 1 : length(sim_t)
    i
    if i <=100
        t_reach = find_earliest_BRS_ind(g_grid, V_R, x(:,i), 101, 1);
        if t_reach>i
            x(:,i+1) = x(:,i);
        else
            dVdt_R = -(eval_u(g_grid,V_R(:,:,t_reach),x(:,i))-...
                eval_u(g_grid,V_R(:,:,t_reach+1),x(:,i)))/dt;

            exrReach.dVdt = dVdt_R;
            exrReach.x = x(:,i);
            exrReach.dRange = dRange;
            exrReach.threshold = [-0.1,0.01];
            [a_r,b_r] = findACS_D(V_R(:,:,i),g_grid,f_dyn,g_dyn,gd_dyn,exrReach);
            [lb_r,ub_r] = findUpLowBound(a_r,b_r,-1,1);
            % exrAvoid.dVdt = dVdt_A;
            exrAvoid.x = x(:,i);
            exrAvoid.dRange = dRange;
            [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrAvoid);
            [lb_a,ub_a] = findUpLowBound(a_a,b_a,-1,1);
            [lb,ub] = CombineBounds([lb_r,lb_a],[ub_r,ub_a]);
            u(i) = (lb+ub)/2 + (rand(1)-0.5)*(ub-lb);
            % u(i) = -1;
            d(i) = 2*(rand-0.5)/5;
            if isnan(u(i))
                break
            end
            [~, x_temp] = ode45(@(t, s) DI_dyn(t, s, u(i),d(i)), [t t+dt], x(:,i));
            x(:,i+1) = x_temp(end,:);
            t = t+dt;
        end
    elseif i > 100 && i<= 200
        t_reach = find_earliest_BRS_ind(g_grid, V_R2, x(:,i), 101, 1)
        if t_reach>i-100
            x(:,i+1) = x(:,i);
        else
            dVdt_R = -(eval_u(g_grid,V_R2(:,:,t_reach),x(:,i))-...
                eval_u(g_grid,V_R2(:,:,t_reach+1),x(:,i)))/dt;

            exrReach.dVdt = dVdt_R;
            exrReach.x = x(:,i);
            exrReach.dRange = dRange;
            exrReach.threshold = [-0.1,0.01];
            [a_r,b_r] = findACS_D(V_R2(:,:,i-100),g_grid,f_dyn,g_dyn,gd_dyn,exrReach);
            [lb_r,ub_r] = findUpLowBound(a_r,b_r,-1,1);
            % exrAvoid.dVdt = dVdt_A;
            exrAvoid.x = x(:,i);
            exrAvoid.dRange = dRange;
            [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrAvoid);
            [lb_a,ub_a] = findUpLowBound(a_a,b_a,-1,1);
            [lb,ub] = CombineBounds([lb_r,lb_a],[ub_r,ub_a]);
            u(i) = (lb+ub)/2 + (rand(1)-0.5)*(ub-lb);

            if isnan(u(i))
                break
            end
            d(i) = 2*(rand-0.5)/5;
            [~, x_temp] = ode45(@(t, s) DI_dyn(t, s, u(i),d(i)), [t t+dt], x(:,i));
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

        [lb_clvf,ub_clvf] = findUpLowBound(a_clvf,b_clvf,-1,1);
        [lb_a,ub_a] = findUpLowBound(a_a,b_a,-1,1);

        [lb,ub] = CombineBounds([lb_clvf,lb_a],[ub_clvf,ub_a]);

        d(i) = 2*(rand-0.5)/5;
        u(i) = (lb+ub)/2 + (rand(1)-0.5)*(ub-lb);
        if isnan(u(i))
                break
        end
        [~, x_temp] = ode45(@(t, s) DI_dyn(t, s, u(i),d(i)), [t t+dt], x(:,i));
        x(:,i+1) = x_temp(end,:);
        t = t+dt;

    end


end

toc


%%

figure
hold on
plot(x(1,:),x(2,:),'k')
visSetIm(g_grid,goal,'m',0);
visSetIm(g_grid,goal2,'m',0);

visSetIm(g_grid,obs,'r',0);
visSetIm(g_grid,V_avoid,'r',0);
visSetIm(g_grid,V_clvf,'k',0.28)
visSetIm(g_grid,V_reach2(:,:,end),'b',0);
visSetIm(g_grid,V_reach(:,:,end),'b',0);
% visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0);

% plot(g_grid.xs{1}(ind_infeasible),g_grid.xs{2}(ind_infeasible),'k.')
% plot(g_grid.xs{1}(ind_nan),g_grid.xs{2}(ind_nan),'r.')


%%
function dydt = DI_dyn(t,y,u,d)
dydt = [y(2)+d;  u];
end
