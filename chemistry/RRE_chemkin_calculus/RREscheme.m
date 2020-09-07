function [tmp]=RREscheme(x,y,omega_i, omega_t, delta_1, delta_t, mu_2, mu_3, varargin)

% varargin holds the color of the plot if needed
col = varargin{1};

% data
alpha_2 = mu_2-delta_1;
alpha_3 = mu_3-delta_t;
delta_2 = delta_t-delta_1;
x = x*100;
y = y*100;

% Plot referential
rep = plot([-20, 30],[0,0], 'k', [0,0],[-20,30], 'k');
hold on

% Plot particle
part = plot(x,y);
part.Marker = '*';
part.MarkerSize = 10;
part.Color = col;

% Plot 1st shock
l_inc = 'Incident shock';
inc = plot([0,30],[0,30*tan(omega_i)],'r');
inc.LineWidth = 2.5;

% Plot expansion fan
i = (alpha_2-alpha_3)/5;
yell = [0.8 0.8 0];
l_exp = 'Expansion fan';
for al = alpha_3:i:alpha_2
    exp = plot([0,-20],[0,20*tan(al)]);
    exp.Color = yell;
end

% Plot interface
l_int = 'Gas interface';
int1 = plot([0,30],[0,0],'c');
int1.LineWidth = 2.5;
% Plot deviated interface
int2 = plot([0,-20],[0,-20*tan(delta_t)],'c');
int2.LineWidth = 2.5;
% Plot transmitted shock
gr = [0 0.5 0];
l_trans = 'Transmitted shock';
trans = plot([0,-20],[0,-20*tan(omega_t)]);
trans.Color = gr;
trans.LineWidth = 2.5;

% Plot trajectory in phase I
l_traj = strcat('Trajectory of the particle in phase I');
traj1 = plot([x,y/tan(omega_i)],[y,y],'--');
traj1.Marker = '*';
traj1.MarkerSize = 10;
traj1.LineWidth = 2.5;
traj1.Color = col;

dist1 = y*sin(mu_2+omega_i-delta_1)/(sin(mu_2)*sin(omega_i));
xprime = y/tan(omega_i)-cos(delta_1)*dist1;
yprime = y-sin(delta_1)*dist1;

traj2 = plot([y/tan(omega_i), xprime],[y, yprime],'--');
traj2.Marker = '*';
traj2.MarkerSize = 10;
traj2.LineWidth = 2.5;
traj2.Color = col;

dist2 = y*sin(mu_3-delta_t+omega_i)/(sin(mu_3-delta_2)*sin(omega_i));
xscd = y/tan(omega_i)-cos(delta_1)*dist2;
yscd = y-sin(delta_1)*dist2;
xtrc = 20+xscd;
ytrc = -xtrc*tan(delta_t);
traj3 = plot([xscd, -20],[yscd,yscd+ytrc],'--');
traj3.LineWidth = 2.5;
traj3.MarkerSize = 10;
traj3.Color = col;

% Plot trajectory in phase II
l_flow = 'Flow motion in the inert phase';
flow1 = plot([30,-8/tan(omega_t)],[-8,-8],'-.m');
flow2 = plot([-8/tan(omega_t),-20],[-8, -8-(20-8/tan(omega_t))*tan(delta_t)], '-.m');
flow1.LineWidth = 2.5;
flow2.LineWidth = 2.5;

% Frame
axis equal
xlim([-20,30])
ylim([-15,15])

% Axes
xlabel("x-axis (cm)", 'FontSize', 18)
ylabel("y-axis (cm)", 'FontSize', 18)
%grid on

% Legends
lgd=legend([inc,exp,int1, trans,traj1,flow1],{l_inc,l_exp,l_int,l_trans, l_traj, l_flow});
legend('boxoff')
lgd.Location = 'northwest';
lgd.FontSize=14;
tmp = true;
end