function [Result] = PrefactorCalculation(PartNum, TrajDiv, R_param, Alpha, eta, Z, Trajectories)
    particle_n = PartNum;
    N_division = TrajDiv;
    z_time = Z;
    for ii = 1:length(Alpha)
        r       = R_param(ii);
        alpha   = Alpha(ii);
        
        trajectory = reshape(Trajectories(ii, :, :), PartNum, TrajDiv);
        %In rder to get the tangent vector first calc the derivative of the curve
        velocity = f_trajectory_diff(N_division, trajectory, z_time);
        
        %using the previous differentiated curve to form normed vectors
        unit_velocity_e = f_velocity_unit_vector(velocity);
        
        %creating a perpendicular vector to e named f
        %Normalization and Nominator are here for error checking purposes
        [f, delta_e, Normalization, Nominator] = f_f_vector(unit_velocity_e);
        
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
        
        %using the rotation matrix above we can create the tau vector basis that
        %goes along the trajectory
        tauspace = f_tau_space(Rot_mtx, Tau01);
        
        %determining the curvature of the trajectory
        C = f_curvature(unit_velocity_e, tauspace, z_time, velocity);

        %B matrix
        B_mtx = f_B_matrix(trajectory, alpha, tauspace, eta);

        %T matrix
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
        [propagator, Trace, Trace2] = f_prefactor( Xi_heun, EigVal, z_time, z_time, alpha, r);
        Result(ii) = f_action(eta, alpha, trajectory, r, z_time, EigVal, propagator);
        disp(['splitting: 1D * sqrt_term * exp(-action) = ', num2str(f_action(eta, alpha, trajectory, r, z_time, EigVal, propagator)) '  ' 'Result = ' num2str(Result(ii))]) 
    end

end

