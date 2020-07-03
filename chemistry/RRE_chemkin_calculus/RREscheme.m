function [tmp]=RREscheme(x,y,omega_i, omega_t, delta_1, delta_t, mu_2, mu_3)

% data
alpha_2 = mu_2-delta_1;
alpha_3 = mu_3-delta_t;
delta_2 = delta_t-delta_1;

% Plot referential
rep = plot([-20, 20],[0,0], 'k', [0,0],[-20,20], 'k');
hold on

% Plot particle
l_part = strcat('Lagrangian particle : x=', num2str(x), 'cm, y=', num2str(y), 'cm');
part = plot(x,y,"*b");

% Plot 1st shock
l_inc = 'Incident shock';
inc = plot([0,20],[0,20*tan(omega_i)],'r');

% Plot expansion fan
i = (alpha_2-alpha_3)/5;
l_exp = 'Expansion fan';
for al = alpha_3:i:alpha_2
    exp = plot([0,-20],[0,20*tan(al)],'y');
end

% Plot interface
l_int = 'Gas interface';
int1 = plot([0,20],[0,0],'c');
% Plot deviated interface
int2 = plot([0,-20],[0,-20*tan(delta_t)],'c');

% Plot transmitted shock
l_trans = 'Transmitted shock';
trans = plot([0,-20],[0,-20*tan(omega_t)],'g');


% Plot trajectory in phase I
l_traj = 'Trajectory of the particle';
traj1 = plot([x,y/tan(omega_i)],[y,y],'-.b');

dist1 = y*sin(mu_2+omega_i-delta_1)/(sin(mu_2)*sin(omega_i));
xprime = y/tan(omega_i)-cos(delta_1)*dist1;
yprime = y-sin(delta_1)*dist1;

traj2 = plot([y/tan(omega_i), xprime],[y, yprime],'-.b');

dist2 = y*sin(mu_3-delta_t+omega_i)/(sin(mu_3-delta_2)*sin(omega_i));
xscd = y/tan(omega_i)-cos(delta_1)*dist2;
yscd = y-sin(delta_1)*dist2;
xtrc = 20+xscd;
ytrc = -xtrc*tan(delta_t);
traj3 = plot([xscd, -20],[yscd,yscd+ytrc],'-.b');

% Plot trajectory in phase II
l_flow = 'Flow motion in the inert phase';
flow1 = plot([20,-y/tan(omega_t)],[-y,-y],'-.m');
flow2 = plot([-y/tan(omega_t),-20],[-y, -y-20*tan(delta_t)], '-.m');

% Frame
axis equal
xlim([-20,20])
ylim([-20,20])

% Axes
xlabel("x-axis")
ylabel("y-axis")

% Legends
legend([part,inc,exp,int1, trans,traj1,flow1],{l_part,l_inc,l_exp,l_int,l_trans, l_traj, l_flow})
legend('boxoff')
legend('Location', 'northwest')
tmp = true;
end