function [tmp]=RRRscheme(x,y,omega_i, omega_t, delta_1, delta_t, omega_r, varargin)

% varargin holds the color of the plot if needed
col = varargin{1};

% data
alpha_r = omega_r+delta_1;
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

% Plot reflected shock
yell = [0.8 0.8 0];
l_refl = 'Reflected shock';
refl = plot([0,-20],[0,20*tan(pi-alpha_r)]);
refl.LineWidth = 2.5;
refl.Color = yell;

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
traj1.Color = col;
traj1.LineWidth = 2.5;

dist1 = y*sin(omega_r-omega_i+delta_1)/(sin(omega_r)*sin(omega_i));
xprime = y/tan(omega_i)-cos(delta_1)*dist1;
yprime = y-sin(delta_1)*dist1;

traj2 = plot([y/tan(omega_i), xprime],[y, yprime],'--');
traj2.Color = col;
traj2.Marker = '*';
traj2.MarkerSize = 10;
traj2.LineWidth = 2.5;

xscd = 20+xprime;
yscd = -xscd*tan(delta_t);
traj3 = plot([xprime, -20],[yprime,yprime+yscd],'--');
traj3.Color = col;
traj3.Marker = '*';
traj3.MarkerSize = 10;
traj3.LineWidth = 2.5;

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
grid on 

% Legends
lgd=legend([inc,refl,int1, trans,traj1,flow1],{l_inc,l_refl,l_int,l_trans, l_traj, l_flow});
legend('boxoff');
lgd.Location = 'northwest';
lgd.FontSize = 14;
tmp = true;
end