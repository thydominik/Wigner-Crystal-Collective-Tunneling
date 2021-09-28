clear all
clc

%loading a trajectory that is previoulsy determined by a MC simulatiton
trajectory_load             = load('200_point_ztime_r_1_3_a_7_eta_20');
% trajectory_load             = load('200_point_ztime_r_2_a_7_eta_20');
% trajectory_load             = load('200_point_ztime_r_1_a_7_eta_20');
% trajectory_load             = load('200_point_ztime_r_1_05_a_7_eta_20');
% trajectory_load             = load('200_point_ztime_r_2_a_5_eta_20');
% trajectory_load             = load('200_point_ztime_r_0_6_a_13_eta_20');
% trajectory_load             = load('200_point_ztime_r_1_3_a_13_eta_20');

%this must be the same eq_pos file that used in the MC, but therefore it's
%not important bc the trajectory's endpoints will be the eq_positions!
equilibrium_positions       = load('eq_pos4_15'); %4th column is the alpha value!
trajectory                  = trajectory_load.position;

%this will select the desired endpoints and alpha value
state                       = 4;           
eq_pos                      = equilibrium_positions.eqpos(:,state);

%it's useful not to set in stone the division and particle number for later
%more-particle cases
[particle_n, N_division]    = size(trajectory);

disp(['particle: ',num2str(particle_n)])
disp(['divisions: ', num2str(N_division)])

%z imaginery time paramter
eps             = 10^-15;                   %this have to be changed manually if it changes in the trajectory code!!!
r               = 1.3;                        %match this with the M.C. simulation
alpha           = eq_pos(4); disp(['Alpha= ', num2str(alpha)])                %from eq_pos file!!!!
eta             = 20;                   %don't try to change this. Fitted value(experimetnts) = 18.813
limits          = 50;       %1.35;          % +/- T
N               = 2000;

%at almost all places the full z time paramterization can be used, but just
%in case I create a reduced z time parameter
z_time_reduced  = linspace(-1 + eps,1 - eps,N_division);
z_time          = linspace(-1, 1, N_division);

figure(1)
clf(figure(1))
hold on
title('3 particle trajectories')
scatter(z_time , trajectory(1,:))
plot(z_time , trajectory(2,:))
plot(z_time , trajectory(3,:))
xlabel('z ')
ylabel('q(z)')
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

%In rder to get the tangent vector first calc the derivative of the curve
velocity = f_trajectory_diff(N_division, trajectory, z_time);

figure(3)
clf(figure(3))
hold on
title('3 velocity curves')
plot(z_time, velocity(1,:))
plot(z_time, velocity(2,:))
plot(z_time, velocity(3,:))
set(gca, 'YScale', 'log')
ylabel('v_0(z) = \chi^\prime (y)')
xlabel('z')
hold off

%using the previous differentiated curve to form normed vectors
unit_velocity_e = f_velocity_unit_vector(velocity);

figure(4)
clf(figure(4))
hold on
title('unit vektor norma check')
ylabel('norm(e)')
xlabel('z')
plot(z_time, unit_velocity_e(1,:).^2 + unit_velocity_e(2,:).^2 + unit_velocity_e(3,:).^2)
set(gca,'Yscale','log')
hold off

figure(5)
clf(figure(5))
hold on
title('3D velocity space')
plot3(velocity(1,:), velocity(2,:), velocity(3,:))
xlabel('v_{q1}')
ylabel('v_{q2}')
zlabel('v_{q3}')
hold off

%creating a perpendicular vector to e named f
%Normalization and Nominator are here for error checking purposes
[f, delta_e, Normalization, Nominator] = f_f_vector(unit_velocity_e);

figure(6)
clf(figure(6))
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

%Error check, how perpendicular the two vectors created above:
for i = 1:N_division
    scalarproduct(i) = f(:,i)'*unit_velocity_e(:,i);
end
figure(7)
clf(figure(7))
hold on
title('scalar product: abs(f * e)')
xlabel('y')
ylabel('f * e')
plot(z_time, abs(scalarproduct))
xlim([z_time(1) 0])
hold off

%determining the angle between consecutive e vectors
delta_phi = zeros(1,N_division);
for i = 1:N_division
    if i == 1
        delta_phi(i) = 2 * asin((norm(delta_e(:,i))/2));
    else
        delta_phi(i) = 2 * asin((norm(delta_e(:,i))/2));
    end
end

figure(8)
clf(figure(8))
hold on
xlabel('y')
ylabel('\Delta \rho(deg°)')
title('\Delta \rho (\tau)')
plot(z_time, rad2deg(delta_phi))
hold off

%rotation matrix, that transforms e_i & f_i to e_i+1 and f_i+1
Rot_mtx = f_rotation_matrix(N_division, unit_velocity_e, f, delta_phi);

%Determining the basis vectors from the equilibrium positions and
%potentials
[Tau01, EigVal] = f_init_tau(trajectory, alpha, velocity, eta);

x0 = trajectory(1,1);
y0 = trajectory(2,1);
z0 = trajectory(3,1);

figure(9)
clf(figure(9))
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

%using the rotation matrix above we can create the tau vector basis that
%goes along the trajectory
tauspace = f_tau_space(Rot_mtx, Tau01);

figure(10)
clf(figure(10))
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

%determining the curvature of the trajectory
C = f_curvature(unit_velocity_e, tauspace, z_time, velocity);

figure(11)
clf(figure(11))
hold on
title('Curvature(z)','FontSize',20)
xlabel('z','FontSize',20)
ylabel('C_{\alpha} (z)','FontSize',20)
lim = 1;
plot(z_time(lim:end-lim), C(1,lim:end-lim),'r')
plot(z_time(lim:end-lim), C(2,lim:end-lim),'k')
plot(z_time(lim:end-lim), C(3,lim:end-lim),'k')
%set(gca,'Yscale','log')
xlim([z_time(1) z_time(end/2)])
%ylim([-0.1 0.1])
hold off


%B matrix
B_mtx = f_B_matrix(trajectory, alpha, tauspace, eta);

figure(19)
clf(figure(19))
hold on
title('B matrix elements')
plot(reshape(B_mtx(1,1,:), [1,200]),'r')
plot(reshape(B_mtx(1,2,:), [1,200]))
plot(reshape(B_mtx(2,1,:), [1,200]))
plot(reshape(B_mtx(2,2,:), [1,200]),'b')
hold off

%T matrix
T_mtx = f_T_matrix(tauspace);

%omega squared part two: first I'll try to do the matrix with the curvatures
%this is already the squared omega matrix hence the '_sq'
omega_sq = f_omega_squared(T_mtx, B_mtx, velocity, C);     



omega_sq_old = omega_sq;

figure(12)
clf(figure(12))
hold on
title('elements of the Omega matrix')
for i = 1:200
    a(i) = omega_sq(1,1,i);
    b(i) = omega_sq(2,1,i);
    c(i) = omega_sq(1,2,i);
    d(i) = omega_sq(2,2,i);
end
plot(z_time, (a(:)),'b')
plot(z_time, (b(:)),'r')
plot(z_time, (c(:)),'ko')
plot(z_time, (d(:)), 'o')
%xlim([-1 0])
hold off

omega_sq = f_fitting(omega_sq, z_time);

figure(13)
clf(figure(13))
hold on
title('elements of the Omega matrix fitted')
for i = 1:200
    aa(i) = omega_sq(1,1,i);
    bb(i) = omega_sq(2,1,i);
    cc(i) = omega_sq(1,2,i);
    dd(i) = omega_sq(2,2,i);
end
plot(z_time, (aa(:)),'b')
plot(z_time, (bb(:)),'r')
plot(z_time, (cc(:)),'ko')
plot(z_time, (dd(:)), 'o')
%xlim([-1 0])
hold off

figure(14)
clf(figure(14))
hold on
title('Fitted (---) and original (o) values of the omega square matrix')
plot(z_time(1:50), (a(1:50)),'bo')
plot(z_time(1:50), (d(1:50)),'ro')
% plot(z_time, (c(:)),'ko')
% plot(z_time, (d(:)), 'o')
plot(z_time(1:50), (aa(1:50)),'b')
plot(z_time(1:50), (dd(1:50)),'r')
% plot(z_time, (cc(:)),'ko')
% plot(z_time, (dd(:)), 'o')
%xlim([-1 0])
hold off

%now only the heun method works, but later on I will be interested in the
%simpler euler method and their difference

[Xi_heun, temp_mtx] = f_diff_equation(z_time, omega_sq, r);

figure(15)
clf(figure(15))
hold on
title('difference between the Xi matrix elements and the omega_0 matrix elements')
for i = 1:198
    temp_a(i) = temp_mtx(1,1,i);
    temp_b(i) = temp_mtx(2,1,i);
    temp_c(i) = temp_mtx(1,2,i);
    temp_d(i) = temp_mtx(2,2,i);
end
plot(z_time(2:end-1), temp_a(:),'b')
plot(z_time(2:end-1), temp_b(:),'r')
plot(z_time(2:end-1), temp_c(:),'ko')
plot(z_time(2:end-1), temp_d(:),'o')
hold off
%%
figure(16)
clf(figure(16))
hold on
title('Xi_{heun} method from the differential equation')
for i = 1:198
    temp_a(i) = Xi_heun(1,1,i);
    temp_b(i) = Xi_heun(2,1,i);
    temp_c(i) = Xi_heun(1,2,i);
    temp_d(i) = Xi_heun(2,2,i);
end
plot(z_time(2:end-1), temp_a(:),'k', 'Displayname','Xi(1,1)')
plot(z_time(2:end-1), temp_b(:),'bo', 'Displayname','Xi(2,1)')
plot(z_time(2:end-1), temp_c(:),'b', 'Displayname','Xi(1,2)')
plot(z_time(2:end-1), temp_d(:),'r', 'Displayname','Xi(2,2)')
xlabel('z')
ylabel('Matrix elements')
hold off
%%
figure(17)
clf(figure(17))
hold on
for i = 1:198
    temp_a(i) = temp_mtx(1,1,i);
    temp_b(i) = temp_mtx(1,2,i);
    temp_c(i) = temp_mtx(2,1,i);
    temp_d(i) = temp_mtx(2,2,i);
end
plot(temp_a(:), 'b')
plot(temp_d(:), 'r')
plot(temp_b(:), 'bo')
plot(temp_c(:), 'ro')
hold off

%%

%constructing the full propagotor with the 1D term
[propagator, Trace, Trace2] = f_prefactor( Xi_heun, EigVal, z_time, z_time, alpha, r);

figure(18)
clf(figure(18))
hold on
title('Trace from T_0 to 0')
plot(z_time(1:end/2), (Trace))
plot(z_time, zeros(1,length(z_time)))
plot(z_time, 1./(1 - z_time.^2))
plot(z_time(1:end/2), 1./(1 - z_time(1:end/2).^2) .* Trace)
%scatter(y_time, zeros(1,length(y_time)))
xlim([z_time(1) 0])
ylim([0 2])
hold off

figure(20)
clf(figure(20))
hold on
title('Trace from T_0 to 0')
plot(z_time(1:end/2), (Trace2),'k','Displayname', 'forward trace')
plot(z_time(1:end/2), (Trace),'r','Displayname', 'backward trace')
plot(z_time(1:end/2), (Trace2 + Trace)/2,'Displayname', 'average trace')
%scatter(y_time, zeros(1,length(y_time)))
xlim([z_time(1) 0])
ylim([0 2])
xlabel('z')
ylabel('Tr(\Omega_0 - \Xi (z))')
hold off

figure(21)
clf(figure(21))
title('determinant of xI')
hold on
for i = 1:200
    deter(i) = det(Xi_heun(:,:,i));
    
end
plot(deter)
hold off

disp(['splitting: 1D * sqrt_term * exp(-action) = ', num2str(f_action(eta, alpha, trajectory, r, z_time, EigVal, propagator))])


