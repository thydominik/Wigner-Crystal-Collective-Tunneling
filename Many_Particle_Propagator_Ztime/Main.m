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

velocity = f_trajectory_diff(N_division, trajectory, z_time);

figure(4)
clf(figure(4))
hold on
title('3 velocity curves')
plot(z_time, velocity(1,:))
plot(z_time, velocity(2,:))
plot(z_time, velocity(3,:))
set(gca, 'YScale', 'log')
ylabel('v_0(z) = \chi^\prime (y)')
xlabel('z')
hold off

%unit velocity

unit_velocity_e = f_velocity_unit_vector(velocity);

figure(5)
clf(figure(5))
hold on
title('unit vektor norma check')
ylabel('norm(e)')
xlabel('z')
plot(z_time, unit_velocity_e(1,:).^2 + unit_velocity_e(2,:).^2 + unit_velocity_e(3,:).^2)
set(gca,'Yscale','log')
hold off

figure(6)
clf(figure(6))
hold on
title('3D velocity space')
plot3(velocity(1,:), velocity(2,:), velocity(3,:))
xlabel('v_{q1}')
ylabel('v_{q2}')
zlabel('v_{q3}')
hold off

%f vector calculation

 [f, delta_e, Normalization, Nominator] = f_f_vector(unit_velocity_e);

figure(7)
clf(figure(7))
hold on
indexing = 2;
title('e & f in v_q space')
plot3(velocity(1,:), velocity(2,:), velocity(3,:),'LineWidth',2)
quiver3(velocity(1,1:indexing:end), velocity(2,1:indexing:end), velocity(3,1:indexing:end), (unit_velocity_e(1,1:indexing:end)+velocity(1,1:indexing:end)), (unit_velocity_e(2,1:indexing:end)+velocity(2,1:indexing:end)), (unit_velocity_e(3,1:indexing:end)+velocity(3,1:indexing:end)),'k');
quiver3(velocity(1,1:indexing:end), velocity(2,1:indexing:end), velocity(3,1:indexing:end), (f(1,1:indexing:end)+velocity(1,1:indexing:end)), (f(2,1:indexing:end)+velocity(2,1:indexing:end)), (f(3,1:indexing:end)+velocity(3,1:indexing:end)),'r');
xlabel('v_{q1}')
ylabel('v_{q2}')
zlabel('v_{q3}')
hold off

for i = 1:N_division
    scalarproduct(i) = f(:,i)'*unit_velocity_e(:,i);
end
figure(8)
clf(figure(8))
hold on
title('scalar product: abs(f * e)')
xlabel('y')
ylabel('f * e')
plot(z_time, abs(scalarproduct))
xlim([z_time(1) 0])
hold off

delta_phi = zeros(1,N_division);
for i = 1:N_division
    if i == 1
        delta_phi(i) = 2 * asin((norm(delta_e(:,i))/2));
    else
        delta_phi(i) = 2 * asin((norm(delta_e(:,i))/2));
    end
end

figure(9)
clf(figure(9))
hold on
xlabel('y')
ylabel('\Delta \rho(deg°)')
title('\Delta \rho (\tau)')
plot(z_time, rad2deg(delta_phi))
hold off

%rotation matrix
%--------------------------------------------------------------------------
Rot_mtx = f_rotation_matrix(N_division, unit_velocity_e, f, delta_phi);
%orthogonal vector basis at -T calculation (\tau)
%Mivel a kezdőpontok segítségével lett meghatározva a pálya, így azok
%végpontjai adják a kezdőpontokat.


[Tau01, EigVal] = f_init_tau(trajectory, alpha, velocity);

x0 = trajectory(1,1);
y0 = trajectory(2,1);
z0 = trajectory(3,1);

figure(10)
clf(figure(10))
hold on
oszto   = 10;
szorzo  = 5*10^7;
plot3(trajectory(1,:) - x0, trajectory(2,:) - y0, trajectory(3,:) - z0,'LineWidth',2)
quiver3(0, 0, 0, (Tau01(1,1))/oszto, (Tau01(2,1))/oszto, (Tau01(3,1))/oszto,'r','LineWidth',2)
quiver3(0, 0, 0, (Tau01(1,2))/oszto, (Tau01(2,2))/oszto, (Tau01(3,2))/oszto,'k','LineWidth',2)
quiver3(0, 0, 0, (Tau01(1,3))/oszto, (Tau01(2,3))/oszto, (Tau01(3,3))/oszto,'k','LineWidth',2)
xlabel('x')
xlim([0 0.05])
ylabel('y')
ylim([0 0.05])
zlabel('z')
zlim([0 0.05])
hold off
%\tau space with roation matrix
%--------------------------------------------------------------------------
tauspace = f_tau_space(Rot_mtx, Tau01);

figure(11)
clf(figure(11))
hold on
s = 10;                  %round(2*N_division/5);
m = 2;
fin = length(trajectory)-s+1;
tau1 = tauspace(:,1,:);
tau2 = tauspace(:,2,:);
tau3 = tauspace(:,3,:);

plot3(trajectory(1,:) - x0, trajectory(2,:) - y0, trajectory(3,:) - z0,'LineWidth',4)
quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau1(1,s:m:fin), tau1(2,s:m:fin), tau1(3,s:m:fin),'r', 'LineWidth',1)
quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau2(1,s:m:fin), tau2(2,s:m:fin), tau2(3,s:m:fin),'k', 'LineWidth',1)
quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau3(1,s:m:fin), tau3(2,s:m:fin), tau3(3,s:m:fin),'k', 'LineWidth',1)
hold off
     
C = f_curvature(unit_velocity_e, tauspace, z_time, velocity);

figure(13)
clf(figure(13))
hold on
title('Curvature(z)','FontSize',20)
xlabel('z','FontSize',20)
ylabel('C_{\alpha} (z)','FontSize',20)
lim = 10;
plot(z_time(lim:end-lim), C(1,lim:end-lim),'r')
plot(z_time(lim:end-lim), C(2,lim:end-lim),'k')
plot(z_time(lim:end-lim), C(3,lim:end-lim),'k')
%set(gca,'Yscale','log')
xlim([z_time(1) z_time(end)])
%ylim([-0.1 0.1])
hold off

%B matrix
%--------------------------------------------------------------------------
B_mtx = f_B_matrix(trajectory, alpha, tauspace);
%T matrix
%--------------------------------------------------------------------------
T_mtx = T_matrix(tauspace);

%omega squared part two: first I'll try to do the matrix with the curvatures
%--------------------------------------------------------------------------
omega_sq = omega_squared(T_mtx, B_mtx, velocity, C);     %this is already the squared omega matrix hence the '_sq'

