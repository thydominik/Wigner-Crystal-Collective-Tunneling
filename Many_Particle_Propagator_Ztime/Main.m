clear all
clc

%we need a trajectory
trajectory_load             = load('200point_ztime');
equilibrium_positions       = load('eq_pos'); %4th column is the alpha value!
trajectory                  = trajectory_load.position;
eq_pos                      = equilibrium_positions.eqpos(:,60);
[particle_n, N_division]    = size(trajectory);

disp(['particle: ',num2str(particle_n)])
disp(['divisions: ', num2str(N_division)])

%z imaginery time paramter:
eps             = 10^-15;                        %this have to be changed manually if it changes in the trajectory code!!!
r               = 3;                        %match this with the M.C. simulation
alpha           = eq_pos(4);                %from eq_pos file!!!!
eta             = 18.813;                   %don't try to change this.
limits          = 50;       %1.35;          % +/- T
N               = 2000;
z_time_reduced  = linspace(-1 + eps,1 - eps,N_division);
z_time          = linspace(-1, 1, N_division);

figure(1)
clf(figure(1))
hold on
title('3 particle trajectories')
scatter(z_time , trajectory(1,:))
plot(z_time , trajectory(2,:))
plot(z_time , trajectory(3,:))
%xlabel('\tau \equiv y = atanh(z)\rho')
%ylabel('q(\tau) \equiv q(y)')
%xlim([-10 10])
hold off

figure(2)
clf(figure(2))
hold on
title('3D coordinate space')
traj1 = trajectory(1,:);
traj2 = trajectory(2,:);
traj3 = trajectory(3,:);
plot3(traj1, traj2, traj3)
xlabel('q_1')
ylabel('q_2')
zlabel('q_3')
hold off

%trajectory derivative


