clear all
clc
disp('-----------------------------------------------------------')
%kell majd egy load trajectory
trajectory_load             = load('traj_100_points');
equilibrium_positions       = load('eq_pos'); %4th column is the alpha value!
trajectory                  = trajectory_load.position;
eq_pos                      = equilibrium_positions.eqpos(:,59);
[particle_n, N_division]    = size(trajectory);

disp(['particle: ',num2str(particle_n)])
disp(['divisions: ', num2str(N_division)])

%z imaginery time paramter:
eps     = 10^-15;                        %this have to be changed manually if it changes in the trajectory code!!!
r       = 3;                        %match this with the M.C. simulation
alpha   = eq_pos(4);                %from eq_pos file!!!!
eta     = 18.813;                   %don't try to change this.
limits  = 50;       %1.35;          % +/- T
N       = 2000;
z_time  = linspace(-1 + eps,1 - eps,N_division);
%insted of this I'll use the not parametrized and dimensionless 'y'
%original imaginery time:
y_time  = atanh(z_time)*r;

figure(1)
clf(figure(1))
hold on
title('3 particle trajectories')
scatter(z_time , trajectory(1,:))
plot(z_time , trajectory(2,:))
plot(z_time , trajectory(3,:))
xlabel('\tau \equiv y = atanh(z)\rho')
ylabel('q(\tau) \equiv q(y)')
%xlim([-10 10])
trajectory2 = trajectory;
y_time2     = y_time;
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
%smooth trajectory vagy illesztés; Függvény: Smoothing 
%--------------------------------------------------------------------------
%% time limit 
Tau0 = init_tau0(eq_pos, alpha);
[new_trajectory,N_division, coefficients] = Smoothing(trajectory, y_time2, limits, N, eq_pos, Tau0);
y_time = linspace(-limits, limits, N_division);         
trajectory = new_trajectory;
%%
x = linspace(-1,1,100);
figure(3)
clf(figure(3))
hold on

scatter(z_time(1:2:end) , trajectory(1,1:2:end),'r')
scatter(z_time(1:2:end) , trajectory(2,1:2:end),'b')
scatter(z_time(1:2:end) , trajectory(3,1:2:end),'k')
xlabel('$z$','interpreter','latex','FontSize',16)
ylabel('$\chi (z)$','interpreter','latex','FontSize',16)
plot(z_time , trajectory2(1,:),'LineWidth',2)
plot(z_time , trajectory2(2,:),'LineWidth',2)
plot(z_time , trajectory2(3,:),'LineWidth',2)

hold off
%%
%trajectroy derivative
%--------------------------------------------------------------------------
velocity = traj_diff(particle_n, N_division, limits, coefficients, eq_pos);

figure(4)
clf(figure(4))
hold on
title('3 velocity curves')
plot(y_time, velocity(1,:))
plot(y_time, velocity(2,:))
plot(y_time, velocity(3,:))
set(gca, 'YScale', 'log')
ylabel('v_0(y) = q^\prime (y)')
xlabel('y')
hold off

%velocity unit vector(e)
%--------------------------------------------------------------------------
unit_velocity_e  = velocity_unit_vector(velocity);

figure(5)
clf(figure(5))
hold on
title('unit vektor norma check')
ylabel('norm(e)')
xlabel('y')
plot(y_time, unit_velocity_e(1,:).^2 + unit_velocity_e(2,:).^2 + unit_velocity_e(3,:).^2)
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
%--------------------------------------------------------------------------
[f, delta_e, Normalization, Nominator] = f_vector(unit_velocity_e);

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
plot(y_time, abs(scalarproduct))
xlim([y_time(1) 0])
hold off

%\phi angle calculation
%--------------------------------------------------------------------------
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
plot(y_time, rad2deg(delta_phi))
hold off
    
%rotation matrix
%--------------------------------------------------------------------------
Rot_mtx = rotation_matrix(delta_phi,unit_velocity_e,f,N_division);

%orthogonal vector basis at -T calculation (\tau)
%Mivel a kezdõpontok segítségével lett meghatározva a pálya, így azok
%végpontjai adják a kezdõpontokat.

[Tau01, Tau02, EigVal] = init_tau(trajectory, alpha, velocity);

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

quiver3(0, 0, 0, (Tau02(1,1))/oszto, (Tau02(2,1))/oszto, (Tau02(3,1))/oszto,'r','LineWidth',2)
quiver3(0, 0, 0, (Tau02(1,2))/oszto, (Tau02(2,2))/oszto, (Tau02(3,2))/oszto,'k','LineWidth',2)
quiver3(0, 0, 0, (Tau02(1,3))/oszto, (Tau02(2,3))/oszto, (Tau02(3,3))/oszto,'k','LineWidth',2)
xlabel('x')
xlim([0 0.05])
ylabel('y')
ylim([0 0.05])
zlabel('z')
zlim([0 0.05])
hold off
%\tau space with roation matrix
%--------------------------------------------------------------------------
tauspace = Tau_space(Rot_mtx, Tau01);
tauspace2= Tau_space(Rot_mtx, Tau02);

figure(11)
clf(figure(11))
hold on
s = 10;                  %round(2*N_division/5);
m = 2;
fin = length(trajectory)-s+1;
tau1 = tauspace2(:,1,:);
tau2 = tauspace2(:,2,:);
tau3 = tauspace2(:,3,:);

plot3(trajectory(1,:) - x0, trajectory(2,:) - y0, trajectory(3,:) - z0,'LineWidth',4)
quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau1(1,s:m:fin), tau1(2,s:m:fin), tau1(3,s:m:fin),'r', 'LineWidth',1)
quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau2(1,s:m:fin), tau2(2,s:m:fin), tau2(3,s:m:fin),'k', 'LineWidth',1)
quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau3(1,s:m:fin), tau3(2,s:m:fin), tau3(3,s:m:fin),'k', 'LineWidth',1)
hold off
     
figure(12)
clf(figure(12))
for i = 1:length(tau1)
    scalaprod12(i) = tau1(:,i)'*tau2(:,i);
    scalaprod13(i) = tau1(:,i)'*tau3(:,i);
    matrixorth(i) = trace(Rot_mtx(:,:,i)*Rot_mtx(:,:,i)');
end
hold on
title('blue: if 0 then its orthogonal mtx; red: tau1*tau2; black: tau1*tau3')
plot(y_time, scalaprod12,'r')
plot(y_time, scalaprod13,'k')
plot(y_time(2:end-1), matrixorth(2:end-1)-3,'b')
xlim([-limits 0])
hold off
    
%curvature calculation
%--------------------------------------------------------------------------
C = curvature_second(unit_velocity_e, tauspace2, y_time, velocity);

figure(13)
clf(figure(13))
hold on
title('Curvature(\tau)','FontSize',20)
xlabel('y \equiv \tau','FontSize',20)
ylabel('C_{\alpha} (\tau)','FontSize',20)
lim = 10;
plot(y_time(lim:end-lim), C(1,lim:end-lim),'r')
plot(y_time(lim:end-lim), C(2,lim:end-lim),'k')
plot(y_time(lim:end-lim), C(3,lim:end-lim),'k')
%set(gca,'Yscale','log')
xlim([y_time(1) y_time(end)])
%ylim([-0.1 0.1])
hold off

%B matrix
%--------------------------------------------------------------------------
B_mtx = B_matrix(trajectory, alpha, tauspace);
%T matrix
%--------------------------------------------------------------------------
T_mtx = T_matrix(tauspace);

%omega squared part two: first I'll try to do the matrix with the curvatures
%--------------------------------------------------------------------------
omega_sq = omega_squared(T_mtx, B_mtx, velocity, C);     %this is already the squared omega matrix hence the '_sq'

%differential equation for captial \xi
%--------------------------------------------------------------------------
[tau_time, xi_mtx] = diff_equation(omega_sq, y_time(1:round(length(y_time)/2)), EigVal);

[Xi_euler, Xi_heun] = diff_equation2(y_time, omega_sq, EigVal);


for i = 1:N_division
    j = 2;
    matrix_difference1(i) = Xi_euler(j,j,i) - Xi_heun(j,j,i);
    matrix_difference2(i) = Xi_euler(j-1,j-1,i) - Xi_heun(j-1,j-1,i);
end


figure(14)
clf(figure(14))
hold on
plot(matrix_difference1)
plot(matrix_difference2)
%set(gca,'Yscale','log')
hold off

%whole prefactor
%--------------------------------------------------------------------------
[propagator1,Trace1] = prefactor( Xi_euler, EigVal, y_time, y_time, alpha);

[propagator2,Trace2] = prefactor( Xi_heun, EigVal, y_time, y_time, alpha);



