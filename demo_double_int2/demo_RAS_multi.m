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

%% Define the system:
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

%%
xs = zeros(prod(g_grid.N, 'all'),2);
xs(:,1) = reshape(g_grid.xs{1}, 1, []);
xs(:,2) = reshape(g_grid.xs{2}, 1, []);

%% Online simulation
dt = 0.01;
sim_t = [0:dt:6];
% x0 = [ -1.1 , 0.4 ];
% x = nan(2,length(sim_t));
% x(:,1) = x0;
t = 0;
ind_verified = [];
n = prod(g_grid.N,'all');
ind_all = [1:prod(g_grid.N,'all')];

xs = zeros(prod(g_grid.N, 'all'),2);
xs(:,1) = reshape(g_grid.xs{1}, 1, []);
xs(:,2) = reshape(g_grid.xs{2}, 1, []);
% xs_bk = xs;


x = cell(2,1);
x{1} = nan(g_grid.N(1),g_grid.N(2),length(sim_t));
x{2} = nan(g_grid.N(1),g_grid.N(2),length(sim_t));

ind_RA_empty = [];
ind_SA_empty = [];
ind_empty = [];

% x = nan(2,length(sim_t));
% u = nan(1,length(sim_t));
% x(:,1) = [g_grid.xs{1}(j),g_grid.xs{2}(j)];
% x(:,1) = x0;
% figure
% hold off
tic
for i = 1 : 100 %length(sim_t)
    i
        t_reach = find_earliest_BRS_ind_ZG(g_grid, V_reach, xs, 101, 1);
        ind_inn = find(t_reach < 101 - i);
        dVdt_R = -(eval_u(g_grid,V_R(:,:,i),xs)-...
            eval_u(g_grid,V_R(:,:,i+1),xs))/dt;
        % dVdt_A = -(eval_u(g_grid,V_A(:,:,i),xs)-...
        %     eval_u(g_grid,V_A(:,:,i+1),xs))/dt;
        exrReach.dVdt = dVdt_R;
        exrReach.x = xs;
        exrReach.dRange = dRange;
        exrReach.threshold = [-0.1,0.01];
        [a_r,b_r] = findACS_D(V_R(:,:,i),g_grid,f_dyn,g_dyn,gd_dyn,exrReach);
        [lb_r,ub_r] = findUpLowBound(a_r,b_r,-1,1);
        % exrAvoid.dVdt = dVdt_A;
        exrAvoid.x = xs;
        exrAvoid.dRange = dRange;
        [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrAvoid);
        [lb_a,ub_a] = findUpLowBound(a_a,b_a,-1,1);
        [lb,ub] = CombineBounds([lb_r,lb_a],[ub_r,ub_a]);
        us = (lb+ub)/2 + (rand([prod(g_grid.N,'all'),1])-0.5).*(ub-lb);
        ds = 0*(rand([prod(g_grid.N,'all'),1])-0.5)/5;
        ind_RA_e = find(isnan(us));
        ind_RA_empty = union(ind_RA_empty,ind_RA_e);
        xs_bk = xs;
        [~,p] = ode45(@(t,p) DI_dyn(t,p,us,ds,n),[t t+dt],xs');
        p = reshape(p,[],2,n);
        xs = squeeze(p(end,:,:))';
        % dx = DI_dyn(t, xs, us,ds);
        % xs = xs+dx*dt;
        xs(ind_inn,:) = xs_bk(ind_inn,:);
        t = t+dt;
        % plot(xs(:,1),xs(:,2),'b.')
        % drawnow
    x{1}(:,:,i) = reshape(xs(:,1),g_grid.N(1),g_grid.N(2));
    x{2}(:,:,i) = reshape(xs(:,2),g_grid.N(1),g_grid.N(2));
end


for i = 1:500
    i
    if i <= 100
        t_reach = find_earliest_BRS_ind_ZG(g_grid, V_reach2, xs, 101, 1);
        ind_inn = find(t_reach < 101 - i );

        dVdt_R = -(eval_u(g_grid,V_R2(:,:,i),xs)-...
            eval_u(g_grid,V_R2(:,:,i+1),xs))/dt;
        exrReach.dVdt = dVdt_R;
        exrReach.x = xs;
        exrReach.dRange = dRange;
        exrReach.threshold = [-0.1,0.01];
        [a_r,b_r] = findACS_D(V_R2(:,:,i),g_grid,f_dyn,g_dyn,gd_dyn,exrReach);
        [lb_r,ub_r] = findUpLowBound(a_r,b_r,-1,1);
        exrAvoid.x = xs;
        exrAvoid.dRange = dRange;

        [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrAvoid);
        [lb_a,ub_a] = findUpLowBound(a_a,b_a,-1,1);
        [lb,ub] = CombineBounds([lb_r,lb_a],[ub_r,ub_a]);
        us = (lb+ub)/2 + (rand([prod(g_grid.N,'all'),1])-0.5).*(ub-lb);
        ds = 2*(rand([prod(g_grid.N,'all'),1])-0.5)/5;
        ind_RA_e = find(isnan(us));
        ind_RA_empty = union(ind_RA_empty,ind_RA_e);
        xs_bk = xs;
        [~,p] = ode45(@(t,p) DI_dyn(t,p,us,ds,n),[t t+dt],xs');
        p = reshape(p,[],2,n);
        xs = squeeze(p(end,:,:))';
        % dx = DI_dyn(t, xs, us,ds);
        % xs = xs+dx*dt;
        xs(ind_inn,:) = xs_bk(ind_inn,:);
        t = t+dt;
    else
        clvf = eval_u(g_grid,V_clvf,xs);
        i_clvf = find(~isnan(clvf));
        if all (clvf(i_clvf)<= 0.5)
            break
        end
        exrCLF.x = xs;
        exrCLF.gamma = 0.1;
        exrCLF.dRange = dRange;
        exrA.gamma = 0;
        exrA.x = xs;
        exrA.dRange = dRange;
        [a_clvf,b_clvf] = findACS_D(V_clvf,g_grid,f_dyn,g_dyn,gd_dyn,exrCLF);
        [a_a,b_a] = findACS_D(V_avoid,g_grid,f_dyn,g_dyn,gd_dyn,exrA);

        [lb_clvf,ub_clvf] = findUpLowBound(a_clvf,b_clvf,-1,1);
        [lb_a,ub_a] = findUpLowBound(a_a,b_a,-1,1);

        [lb,ub] = CombineBounds([lb_clvf,lb_a],[ub_clvf,ub_a]);
    

        us = (lb+ub)/2 + (rand([prod(g_grid.N,'all'),1])-0.5).*(ub-lb);
        ds = 2*(rand([prod(g_grid.N,'all'),1])-0.5)/5;
        ind_SA_e = find(isnan(us));
        [~,p] = ode45(@(t,p) DI_dyn(t,p,us,ds,n),[t t+dt],xs');
        p = reshape(p,[],2,n);
        xs = squeeze(p(end,:,:))';
        % dx = DI_dyn(t, xs, us,ds);
        % xs = xs+dx*dt;
        t = t+dt;
        ind_out = find(xs(:,1)>2 | xs(:,1)<-2 | xs(:,2)>2 | xs(:,2)<-2);
        xs(ind_out,1) = nan;
        xs(ind_out,2) = nan;
       
        ind_SA_empty = union(ind_SA_empty,union(ind_out,ind_SA_e));
    end
    

    x{1}(:,:,i+100) = reshape(xs(:,1),g_grid.N(1),g_grid.N(2));
    x{2}(:,:,i+100) = reshape(xs(:,2),g_grid.N(1),g_grid.N(2));
end
toc
ind_empty = intersect(ind_RA_empty,ind_SA_empty);

%%
ind_verified = setdiff(ind_all,ind_empty);
i = 140;
j =250;

xp = reshape(x{1}(i,j,:),1,[]);
yp = reshape(x{2}(i,j,:),1,[]);

figure
hold on
% plot(g_grid.xs{1}(ind), g_grid.xs{2}(ind),'r.')
plot(g_grid.xs{1}(ind_verified), g_grid.xs{2}(ind_verified),'b.');
visSetIm(g_grid,goal,'m',0);
visSetIm(g_grid,obs,'r',0);
visSetIm(g_grid,V_avoid,'r',0);
visSetIm(g_grid,V_clvf,'k',0.28);
visSetIm(g_grid,V_reach2(:,:,101),'r',0);

visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0);
plot(xp,yp,'k*')

% plot(g_grid.xs{1}(ind_infeasible),g_grid.xs{2}(ind_infeasible),'k.')
% plot(g_grid.xs{1}(ind_nan),g_grid.xs{2}(ind_nan),'r.')

%%
% MakeVideo = 0;
% 
% myfig = figure()
% hold on
% set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.8]);
% hold on
% h_trajV = animatedline('Marker','o');
% visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0);
% visSetIm(g_grid,V_clvf,'r',0.5)
% visSetIm(g_grid,obs,'k',0)
% visSetIm(g_grid,goal,'g',0)
% xlabel('x1')
% ylabel('x2')
% % zlabel('value')
% % title('value function')
% if MakeVideo==1
%     v = VideoWriter('RA_correct_traj','MPEG-4');
%     v.FrameRate = 30;
%     open(v);
% 
% end
% 
% for i = 1:400
% 
% 
%     addpoints(h_trajV,xp(i),yp(i));
%     drawnow
%     frame = getframe(gcf);
%     if MakeVideo == 1
%         writeVideo(v,frame);
%     end
% end 
% 
% if MakeVideo == 1
%     close(v);
% end
% 

%%
function dydt = DI_dyn(t,y,u,d,n)
y = reshape(y,[],n);
d = d';
u = u';

dydt = [y(2,:)+d; u];

dydt = dydt(:);
end
