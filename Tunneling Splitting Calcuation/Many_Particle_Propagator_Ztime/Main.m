clc
clear all

disp('Instanton Prefactor and Tunneling splitting calculation')

% this data variable will hold all results from this calculation
data = zeros(15,2);

%This loops goes through the previously calculated trajectories:
for ind = 1:5:71
    %Clearing a variable in each cycle: This will hold some of the arc
    %length
    Sf = [];  

    % Pulling in the trajectory
    name = ind;
    nameSTR = ['P_' num2str(name)];

    %loading a trajectory that is previoulsy determined by a MC simulatiton
    trajectory_load              = load(nameSTR);%load('200_point_ztime_r_1_3_a_8_eta_20');

    %this must be the same eq_pos file that used in the MC, but therefore it's
    %not important bc the trajectory's endpoints will be the eq_positions!
    equilibrium_positions       = load('EqPos_eta20_alpha_5_20'); %4th column is the alpha value!
    trajectory                  = trajectory_load.position;

    %this will select the desired endpoints and alpha value
    state                       = ind;  
    eq_pos                      = equilibrium_positions.eqpos(:, state);

    %it's useful not to set in stone the division and particle number for later
    %more-particle cases
    [particle_n, N_division]    = size(trajectory);

    disp(['Particle # = ',num2str(particle_n)])
    disp(['Trajectory points = ', num2str(N_division)])

    %This needs to be carried over from the trajectory calculation
    R = linspace(1.3, 0.3, 151);
    
    % z imaginery time paramter
    eps             = 10^-15;   %this have to be changed manually if it changes in the trajectory code!!!
    r               = R(state); %match this with the M.C. simulation
    alpha           = eq_pos(4); disp(['Alpha= ', num2str(alpha)])  %from eq_pos file!!!!
    eta             = 20;       %don't try to change this. Fitted value(experimetnts) = 18.813
    limits          = 50;       %1.35;  % +/- T

    %at almost all places the full z time paramterization can be used, but just
    %in case I create a reduced z time parameter
    z_time_reduced  = linspace(-1 + eps,1 - eps,N_division);
    z_time          = linspace(-1, 1, N_division);

    % Using the trajectories, creating the arc length param and the arc
    % length parametrized V(S) potential.
    [chiS, S, VS]   = f_arclength(trajectory, alpha, eta, z_time);
    
    % Obviously initially S \in [0, 2*S_0], but I want S(\tau = 0 ) = 0 so
    % S = [-S_0, S_0]
    S = S - max(S)/2;

    % New Part for the NEW SCHRÖDINGER ED:----------------------------------------
    %dS              = 0.005;     %The division we wnat the arc length
    %Sq              = min(S):dS:max(S);
    Sq              = linspace(min(S), max(S), 200);
    dS              = Sq(2) - Sq(1);
    LSq             = length(Sq);
    VS_interpolate  = interp1(S, VS, Sq, 'Spline');

    figure(2)
    clf(figure(2))
    hold on
    plot(S, VS - min(VS), '.')
    plot(Sq, VS_interpolate - min(VS_interpolate), '.')
    hold off

    % Fitting the V(S) potential with a quartic potential
    [gof, fc] = f_fitting_VS_2(S(1:20), VS(1:20));
    
    figure(6)
    clf(figure(6))
    hold on
    xx = linspace(-3, 3, 1000);
    %plot(xx, fc.a + fc.b * (xx.^2 - fc.c))
    plot(xx, fc.a + fc.b * (xx + min(S)).^2)
    plot(S, VS)
    hold off

    % NEW PART: GLUING TOGGETHER OF THE V(S) POTENTIAL:
    ex = 0.2;      % Variable that tells how much further we want to go from min(S)
    for i = 1:(ex * LSq)
        Sf(i) = (-i) * dS - max(Sq);
    end
    Sf = [flip(Sf) Sq];
    for i = 1:(ex * LSq)
        Sf(end + 1) = i * dS + max(Sq);
    end
    
    %Construct the whole V(S) potential
    %mini = fc.c;
    %mini = min(S)^2;
    %VSf = [(fc.b * (Sf(1:LSq*ex).^2 - mini).^2) (VS_interpolate - min(VS_interpolate)) (fc.b * (Sf(end - LSq*ex + 1:end).^2 - mini).^2)];
    
    figure(3)
    clf(figure(3))
    hold on
    %plot(Sf, VSf, '.-')
    plot(S, VS - min(VS), 'o-')
    hold off

    % Calculating the frequencies from the fits
    %omegaS_sq       = 4 * fc.b * (3 * max(S)^2 - fc.c);
    omegaS_sq_S     = 4 * fc.b * (3 * max(S)^2 - max(S)^2);
    omegaS_sq = omegaS_sq_S;

    disp('Classical frequancy squared from arc length param from 3 parameter fit: ')
    %disp([num2str(omegaS_sq) ' & ' num2str(sqrt(omegaS_sq))])
    disp('Classical frequancy squared from arc length param: ')
    disp([num2str(omegaS_sq_S) ' & ' num2str(sqrt(omegaS_sq_S)) ])
    disp('Difference between the two freqs : ')
    %disp(num2str(abs(sqrt(omegaS_sq_S) - sqrt(omegaS_sq))))

    % Using the fitted V(S) potential I do an ED in this effectively 1
    % dimensional problem
    %[Spectra]   = Schrodinger_VS(VS, S, fc.a, fc.b, fc.c);
    %[Spectra]   = Schrodinger_VS(VS, S, fc.a, fc.b, fc.c);
    %[Psi2, Spectra2]  = Schrodinger_VSF(Sf, VSf, fc.a);

    %Getting the tunneling splitting from this image:
    %E = Spectra;
%     E2 = Spectra2;
%     dE1(ind) = E(2) - E(1);
%     dE2(ind) = E2(2) - E2(1);

    %In order to get the tangent vector first calc the derivative of the curve
    velocity = f_trajectory_diff(N_division, trajectory, z_time);

    %using the previous differentiated curve to form normed vectors
    unit_velocity_e = f_velocity_unit_vector(velocity);

    %creating a perpendicular vector to e named f
    %Normalization and Nominator are here for error checking purposes
    [f, delta_e, Normalization, Nominator] = f_f_vector(unit_velocity_e);

    %Error check, how perpendicular the two vectors created above:
    for i = 1:N_division
        scalarproduct(i) = f(:,i)' * unit_velocity_e(:,i);
    end

    %determining the angle between consecutive e vectors
    delta_phi = zeros(1,N_division);
    for i = 1:N_division
        if i == 1
            delta_phi(i) = 2 * asin((norm(delta_e(:,i))/2));
        else
            delta_phi(i) = 2 * asin((norm(delta_e(:,i))/2));
        end
    end

    %rotation matrix, that transforms e_i & f_i to e_i+1 and f_i+1
    Rot_mtx = f_rotation_matrix(N_division, unit_velocity_e, f, delta_phi);

    %Determining the basis vectors from the equilibrium positions and
    %potentials
    [Tau01, EigVal] = f_init_tau(trajectory, alpha, velocity, eta);

    x0 = trajectory(1,1);
    y0 = trajectory(2,1);
    z0 = trajectory(3,1);

    %using the rotation matrix above we can create the tau vector basis that
    %goes along the trajectory
    tauspace = f_tau_space(Rot_mtx, Tau01);
    
    %determining the curvature of the trajectory
    C = f_curvature(unit_velocity_e, tauspace, z_time, velocity);
    
    %B matrix -- See notes what this is
    B_mtx = f_B_matrix(trajectory, alpha, tauspace, eta);

    %T matrix -- See notes on what this is
    T_mtx = f_T_matrix(tauspace);

    %omega squared part two: first I'll try to do the matrix with the curvatures
    %this is already the squared omega matrix hence the '_sq'
    omega_sq = f_omega_squared(T_mtx, B_mtx, velocity, C);     

    omega_sq_old = omega_sq;

    omega_sq = f_fitting(omega_sq, z_time);

    %now only the heun method works, but later on I will be interested in the
    %simpler euler method and their difference
    [Xi_heun, temp_mtx] = f_diff_equation(z_time, omega_sq, r);

    %constructing the full propagotor with the 1D term
    [propagator, Trace, Trace2] = f_prefactor(Xi_heun, EigVal, z_time, z_time, alpha, r);


    [splittt, action, ODPM1, ODPL, ODPI, ODPM2, ODPM3] = f_action(eta, alpha, trajectory, r, z_time, EigVal, propagator, S, VS, 1);
    disp(['splitting: 1D * sqrt_term * exp(-action) = ', num2str(splittt)])
    data(ind, 1) = alpha;       %positive apha values
    data(ind, 2) = splittt;     %N-1 dimensional part
    data(ind, 3) = action;      %actio calculated from trajectories
    data(ind, 4) = ODPM1;       %One Dimensional Part 1D Milnikov (without integral)
    data(ind, 5) = ODPL;        %Landau prefactor
    data(ind, 6) = ODPI;        %Instanton
    data(ind, 7) = ODPM2;       %Milnikov with integral and original omega
    data(ind, 8) = ODPM3;       %Milnikov with integral and 'corrected' omega
    disp('---------------------------------000-------------------------------------')
end

    delta                       = [];
    delta(:, 1)                 = data(:, 1);
    delta(:, 2)                 = data(:,2);
    data( ~any(data,2), : )     = [];
    delta( ~any(delta,2), :)    = [];

save('WorkSpace')
    
save('data')

figure(1)
clf(figure(1))
hold on
plot(abs(data(:, 1)), data(:,2),'.-', 'DisplayName', 'N-1 Dimensional part')
%plot(abs(data(:, 1)), nonzeros(EEE4), '.-')
%set(gca, 'Yscale', 'log')
xlim([5 12])
legend
grid on
xlabel('|\alpha|')
hold off

figure(2)
clf(figure(2))
hold on
% plot(abs(data(:,1)), data(:,4),'.-', 'DisplayName', 'W/O integral')
% plot(abs(data(:,1)), data(:,7),'.-', 'DisplayName', 'original omega')
% plot(abs(data(:,1)), data(:,8),'.-', 'DisplayName', 'V(S) Omega')
plot(abs(data(:,1)), data(:,4)./data(:,7),'.-', 'DisplayName', 'R')
%plot(abs(data(:,1)), data(:,4)./data(:,8),'.-', 'DisplayName', 'V(S) omega')
%set(gca, 'Yscale', 'log')
xlim([5 12])
title('Ratio of the 1D parts')
grid on
xlabel('|\alpha|')

hold off

PascuD = load('Delta_E_DMRG_Norb_8_eta_20.00.mat');
PascuD2 = load('E_Schrodinger_3e_eta_20.00_N_100_beta_0.300.dat');
M_data = load('E_Schrodinger_3e_eta_20.00_beta_0.01_N_100.dat');
% Putting a 1/sqrt(pi) factor in the splitting calculation pushed the
% Instanton curves down to the ED pretty well
%%

figure(5)
clf(figure(5))
hold on
title('ED from V(S) & Milnikov 1D part')
plot(-M_data(:, 1), M_data(:, 3) - M_data(:, 2), '*-', 'DisplayName', 'ED')
%plot(abs(data(:,1)), data(:, 4) .* data(:, 2), 'o-', 'DisplayName', '1D Milnikov')
%plot(abs(data(:,1)), data(:, 4) .* data(:, 2), '.-', 'DisplayName', '1D Milnikov w/o Int')
plot(abs(data(:,1)), data(:, 7) .*data(:, 2), 'o-', 'DisplayName', '1D Milnikov w/ Int')
%plot(abs(data(:,1)), data(:, 8) .* data(:, 2), 'o-', 'DisplayName', '1D Milnikov w/ Int and omega from arc length param')
plot(abs(data(:, 1)), nonzeros(dE1) .* data(:, 2), 's-', 'DisplayName', '\DeltaE_{1, 2} from V(S)')
plot(abs(data(:, 1)), nonzeros(dE2) .* data(:, 2), '^-', 'DisplayName', '\DeltaE_{1, 2}')
%plot(abs(data(:,1)), data(:,5), '.-', 'LineWidth', 2, 'DisplayName', 'Landau prefactor w/o N-1 dim part')
%plot(abs(data(:,1)), data(:,6) .* data(:,2), 'x-', 'LineWidth', 1, 'DisplayName', 'Coleman pref.')
%plot(abs(data(:,1)), data(:,7) .* data(:,2), '.-', 'LineWidth', 2, 'DisplayName', 'Milnikov')
%plot(abs(PascuD.alpha_list), PascuD.delta_E_DMRG,'.-', 'DisplayName', 'Pascu DMRG')
%plot(abs(PascuD2(:,1)), abs(PascuD2(:, 2) - PascuD2(:, 10)),'.-', 'DisplayName', 'Pascu Schrödinger')
%set(gca, 'Yscale', 'log')
%xlim([5 12])
legend
grid on
xlabel('|\alpha|')
ylabel('\Delta')
%ylim([10^-5 10])
hold off

%%
figure(4)
clf(figure(4))
hold on
plot(abs(data(:,1))-5, data(:,2),'.-', 'DisplayName', 'Instanton')
%plot(dE(:, 1), dE(:, 2), '.-', 'DisplayName', 'Schrödinger & DMRG')
plot(splitt(:, 1), splitt(:, 2), '.-', 'DisplayName', 'Instanton')
set(gca, 'Yscale', 'log')
xlim([0 10])
legend
grid on
xlabel('|\alpha|')
ylabel('\Delta')
hold off
%%
%     figure(20)
%     clf(figure(20))
%     hold on
%     plot(S/max(S), (VS - min(VS))/max((VS - min(VS))), '.-')
%     legend
%     hold off

%     figure(11)
%     %clf(figure(1))
%     hold on
% 
% %     VS = VS - min(VS);
% %     VS = VS./max(VS);
%     plot(S/max(S), (VS - min(VS))/max((VS - min(VS))), '.-')
%     xlabel('S - symm.')
%     ylabel('V(S)')
%     hold off
    
%     figure(2)
%     clf(figure(2))
%     hold on
%     E = diag(Spectra);
%     plot(E)
%     disp(E(3) - E(1))
%     EEE1(ind) = E(2) - E(1);
%     EEE2(ind) = E(3) - E(1);
%     EEE3(ind) = E(3) - E(2);
%     EEE4(ind) = E(4) - E(1);
%     hold off
    
    % figure(1)
    % clf(figure(1))
    % hold on
    % title('3 particle trajectories')
    % scatter(z_time , trajectory(1,:))
    % plot(z_time , trajectory(2,:))
    % plot(z_time , trajectory(3,:))
    % xlabel('z ')
    % ylabel('q(z)')
    % %xlim([-10 10])
    % hold off

    % figure(2)
    % clf(figure(2))
    % hold on
    % title('3D coordinate space')
    % traj1 = trajectory(1,:);
    % traj2 = trajectory(2,:);
    % traj3 = trajectory(3,:);
    % plot3(traj1, traj2, traj3)
    % xlabel('q_1')
    % ylabel('q_2')
    % zlabel('q_3')
    % hold off

    % figure(3)
    % clf(figure(3))
    % hold on
    % title('3 velocity curves')
    % plot(z_time, velocity(1,:))
    % plot(z_time, velocity(2,:))
    % plot(z_time, velocity(3,:))
    % set(gca, 'YScale', 'log')
    % ylabel('v_0(z) = \chi^\prime (y)')
    % xlabel('z')
    % hold off

    % figure(4)
    % clf(figure(4))
    % hold on
    % title('unit vektor norma check')
    % ylabel('norm(e)')
    % xlabel('z')
    % plot(z_time, unit_velocity_e(1,:).^2 + unit_velocity_e(2,:).^2 + unit_velocity_e(3,:).^2)
    % set(gca,'Yscale','log')
    % hold off

    % figure(5)
    % clf(figure(5))
    % hold on
    % title('3D velocity space')
    % plot3(velocity(1,:), velocity(2,:), velocity(3,:))
    % xlabel('v_{q1}')
    % ylabel('v_{q2}')
    % zlabel('v_{q3}')
    % hold off

    

    % figure(6)
    % clf(figure(6))
    % hold on
    % indexing = 2;
    % title('e & f in v_q space')
    % plot3(velocity(1,:), velocity(2,:), velocity(3,:),'LineWidth',2)
    % quiver3(velocity(1,1:indexing:end), velocity(2,1:indexing:end), velocity(3,1:indexing:end), (unit_velocity_e(1,1:indexing:end)+velocity(1,1:indexing:end)), (unit_velocity_e(2,1:indexing:end)+velocity(2,1:indexing:end)), (unit_velocity_e(3,1:indexing:end)+velocity(3,1:indexing:end)),'k');
    % quiver3(velocity(1,1:indexing:end), velocity(2,1:indexing:end), velocity(3,1:indexing:end), (f(1,1:indexing:end)+velocity(1,1:indexing:end)), (f(2,1:indexing:end)+velocity(2,1:indexing:end)), (f(3,1:indexing:end)+velocity(3,1:indexing:end)),'r');
    % xlabel('v_{q1}')
    % ylabel('v_{q2}')
    % zlabel('v_{q3}')
    % hold off

    
    % figure(7)
    % clf(figure(7))
    % hold on
    % title('scalar product: abs(f * e)')
    % xlabel('y')
    % ylabel('f * e')
    % plot(z_time, abs(scalarproduct))
    % xlim([z_time(1) 0])
    % hold off

    

    % figure(8)
    % clf(figure(8))
    % hold on
    % xlabel('y')
    % ylabel('\Delta \rho(deg°)')
    % title('\Delta \rho (\tau)')
    % plot(z_time, rad2deg(delta_phi))
    % hold off

        % figure(9)
    % clf(figure(9))
    % hold on
    % oszto   = 10;
    % szorzo  = 5*10^7;
    % plot3(trajectory(1,:) - x0, trajectory(2,:) - y0, trajectory(3,:) - z0,'LineWidth',2)
    % quiver3(0, 0, 0, (Tau01(1,1))/oszto, (Tau01(2,1))/oszto, (Tau01(3,1))/oszto,'r','LineWidth',2)
    % quiver3(0, 0, 0, (Tau01(1,2))/oszto, (Tau01(2,2))/oszto, (Tau01(3,2))/oszto,'k','LineWidth',2)
    % quiver3(0, 0, 0, (Tau01(1,3))/oszto, (Tau01(2,3))/oszto, (Tau01(3,3))/oszto,'k','LineWidth',2)
    % xlabel('x')
    % xlim([0 0.05])
    % ylabel('y')
    % ylim([0 0.05])
    % zlabel('z')
    % zlim([0 0.05])
    % hold off

    

    % figure(10)
    % clf(figure(10))
    % hold on
    % s = 10;                  %round(2*N_division/5);
    % m = 2;
    % fin = length(trajectory)-s+1;
    % tau1 = tauspace(:,1,:);
    % tau2 = tauspace(:,2,:);
    % tau3 = tauspace(:,3,:);
    % 
    % plot3(trajectory(1,:) - x0, trajectory(2,:) - y0, trajectory(3,:) - z0,'LineWidth',4)
    % quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau1(1,s:m:fin), tau1(2,s:m:fin), tau1(3,s:m:fin),'r', 'LineWidth',1)
    % quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau2(1,s:m:fin), tau2(2,s:m:fin), tau2(3,s:m:fin),'k', 'LineWidth',1)
    % quiver3(trajectory(1,s:m:fin) - x0, trajectory(2,s:m:fin) - y0, trajectory(3,s:m:fin) - z0,tau3(1,s:m:fin), tau3(2,s:m:fin), tau3(3,s:m:fin),'k', 'LineWidth',1)
    % hold off

    

    % figure(11)
    % clf(figure(11))
    % hold on
    % title('Curvature(z)','FontSize',20)
    % xlabel('z','FontSize',20)
    % ylabel('C_{\alpha} (z)','FontSize',20)
    % lim = 1;
    % plot(z_time(lim:end-lim), C(1,lim:end-lim),'r')
    % plot(z_time(lim:end-lim), C(2,lim:end-lim),'k')
    % plot(z_time(lim:end-lim), C(3,lim:end-lim),'k')
    % %set(gca,'Yscale','log')
    % xlim([z_time(1) z_time(end/2)])
    % %ylim([-0.1 0.1])
    % hold off


    

    % figure(19)
    % clf(figure(19))
    % hold on
    % title('B matrix elements')
    % plot(reshape(B_mtx(1,1,:), [1,200]),'r')
    % plot(reshape(B_mtx(1,2,:), [1,200]))
    % plot(reshape(B_mtx(2,1,:), [1,200]))
    % plot(reshape(B_mtx(2,2,:), [1,200]),'b')
    % hold off

    

    

    % figure(12)
    % clf(figure(12))
    % hold on
    % title('elements of the Omega matrix')
    % for i = 1:200
    %     a(i) = omega_sq(1,1,i);
    %     b(i) = omega_sq(2,1,i);
    %     c(i) = omega_sq(1,2,i);
    %     d(i) = omega_sq(2,2,i);
    % end
    % plot(z_time, (a(:)),'b')
    % plot(z_time, (b(:)),'r')
    % plot(z_time, (c(:)),'ko')
    % plot(z_time, (d(:)), 'o')
    % %xlim([-1 0])
    % hold off

    

    % figure(13)
    % clf(figure(13))
    % hold on
    % title('elements of the Omega matrix fitted')
    % for i = 1:200
    %     aa(i) = omega_sq(1,1,i);
    %     bb(i) = omega_sq(2,1,i);
    %     cc(i) = omega_sq(1,2,i);
    %     dd(i) = omega_sq(2,2,i);
    % end
    % plot(z_time, (aa(:)),'b')
    % plot(z_time, (bb(:)),'r')
    % plot(z_time, (cc(:)),'ko')
    % plot(z_time, (dd(:)), 'o')
    % %xlim([-1 0])
    % hold off

    % figure(14)
    % clf(figure(14))
    % hold on
    % title('Fitted (---) and original (o) values of the omega square matrix')
    % plot(z_time(1:50), (a(1:50)),'bo')
    % plot(z_time(1:50), (d(1:50)),'ro')
    % % plot(z_time, (c(:)),'ko')
    % % plot(z_time, (d(:)), 'o')
    % plot(z_time(1:50), (aa(1:50)),'b')
    % plot(z_time(1:50), (dd(1:50)),'r')
    % % plot(z_time, (cc(:)),'ko')
    % % plot(z_time, (dd(:)), 'o')
    % %xlim([-1 0])
    % hold off

    

    % figure(15)
    % clf(figure(15))
    % hold on
    % title('difference between the Xi matrix elements and the omega_0 matrix elements')
    % for i = 1:198
    %     temp_a(i) = temp_mtx(1,1,i);
    %     temp_b(i) = temp_mtx(2,1,i);
    %     temp_c(i) = temp_mtx(1,2,i);
    %     temp_d(i) = temp_mtx(2,2,i);
    % end
    % plot(z_time(2:end-1), temp_a(:),'b')
    % plot(z_time(2:end-1), temp_b(:),'r')
    % plot(z_time(2:end-1), temp_c(:),'ko')
    % plot(z_time(2:end-1), temp_d(:),'o')
    % hold off

    % figure(16)
    % clf(figure(16))
    % hold on
    % title('Xi_{heun} method from the differential equation')
    % for i = 1:198
    %     temp_a(i) = Xi_heun(1,1,i);
    %     temp_b(i) = Xi_heun(2,1,i);
    %     temp_c(i) = Xi_heun(1,2,i);
    %     temp_d(i) = Xi_heun(2,2,i);
    % end
    % plot(z_time(2:end-1), temp_a(:),'k', 'Displayname','Xi(1,1)')
    % plot(z_time(2:end-1), temp_b(:),'bo', 'Displayname','Xi(2,1)')
    % plot(z_time(2:end-1), temp_c(:),'b', 'Displayname','Xi(1,2)')
    % plot(z_time(2:end-1), temp_d(:),'r', 'Displayname','Xi(2,2)')
    % xlabel('z')
    % ylabel('Matrix elements')
    % hold off

    % figure(17)
    % clf(figure(17))
    % hold on
    % for i = 1:198
    %     temp_a(i) = temp_mtx(1,1,i);
    %     temp_b(i) = temp_mtx(1,2,i);
    %     temp_c(i) = temp_mtx(2,1,i);
    %     temp_d(i) = temp_mtx(2,2,i);
    % end
    % plot(temp_a(:), 'b')
    % plot(temp_d(:), 'r')
    % plot(temp_b(:), 'bo')
    % plot(temp_c(:), 'ro')
    % hold off



    

    % figure(18)
    % clf(figure(18))
    % hold on
    % title('Trace from T_0 to 0')
    % plot(z_time(1:end/2), (Trace))
    % plot(z_time, zeros(1,length(z_time)))
    % plot(z_time, 1./(1 - z_time.^2))
    % plot(z_time(1:end/2), 1./(1 - z_time(1:end/2).^2) .* Trace)
    % %scatter(y_time, zeros(1,length(y_time)))
    % xlim([z_time(1) 0])
    % ylim([0 2])
    % hold off

    % figure(20)
    % clf(figure(20))
    % hold on
    % title('Trace from T_0 to 0')
    % plot(z_time(1:end/2), (Trace2),'k','Displayname', 'forward trace')
    % plot(z_time(1:end/2), (Trace),'r','Displayname', 'backward trace')
    % plot(z_time(1:end/2), (Trace2 + Trace)/2,'Displayname', 'average trace')
    % %scatter(y_time, zeros(1,length(y_time)))
    % xlim([z_time(1) 0])
    % ylim([0 2])
    % xlabel('z')
    % ylabel('Tr(\Omega_0 - \Xi (z))')
    % hold off

    % figure(21)
    % clf(figure(21))
    % title('determinant of xI')
    % hold on
    % for i = 1:200
    %     deter(i) = det(Xi_heun(:,:,i));
    %     
    % end
    % plot(deter)
    % hold off