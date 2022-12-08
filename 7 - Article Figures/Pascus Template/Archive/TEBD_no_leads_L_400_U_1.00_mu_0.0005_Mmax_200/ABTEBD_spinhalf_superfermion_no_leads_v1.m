%% ABelian imTEBD to calcute the  Lindblad evolution of a spinless Hamiltonian
%% The parameters are read from the input_file_spinless
% Clear all variables before starting.

clear all;
close all;

%% load the path for the repository
%local folder where distribution is kept.
%localpwd = '/Users/mocap/Documents/research/Numerical_methods/ABTensor_P';
localpwd = '/home/mocap/Documents/research/ABTensor_P';

%add to the path all the subfolders. Each subfolder may contain a given
%matlab function required for the run.

addpath(genpath(localpwd));

%% Input file. I contains the values of different parameters that can be fixed.
% Parameters can be 'numeric', 'string' or 'cell'

input_file='input_file_spinhalf_superfermion';


%% Parameters section. Read the the parameters from the input file.

%chain length;
L = read_variable('L','numeric', input_file);

%On-site Hubbard term;
U = read_variable('U','numeric', input_file);

%chemical potential
mu = read_variable('mu','numeric', input_file);

%Maximum bond dimension
Mmax_TEBD = read_variable('Mmax', 'numeric', input_file);

%Maximum number of steps performed in the TEBD
iteration_no = read_variable('iteration_no', 'numeric', input_file);

% time step
dt = read_variable('dt', 'numeric', input_file);

final_time = dt*(iteration_no-1);


%to get the NESS we fix a relative large final time.
final_time = round(0.8*L);
iteration_no = floor(final_time/dt)

%measuremet time step
dtm = read_variable('dtm', 'numeric', input_file);

%dtm = L;

ratio = floor(dtm/dt)
iteration_m = floor(iteration_no/ratio)


eps = read_variable('eps', 'numeric', input_file);


%%  Symmetries section.

U1sym = generate_U1_Symmetry(1:21,1:21,1:41);


% read the name of the symmetries
symmetry_names = read_variable('symmetry_name', 'cell', input_file);

% get the number of symmetries.
no_of_symmetries = numel(symmetry_names);
symmetries = cell(1, no_of_symmetries);

% convention used for the creation/annihilation super fermions

convention = read_variable('convention', 'string', input_file);


%All symmetries are of U1 type and are arranged in a cell array.
for symmetry =1:no_of_symmetries
    symmetries{symmetry} = U1sym;
end


%% Construct the quantum numbers for the system.
Q =0;
twoSz = 0;

%% Generating the goal representation based on the twoSz and Q
% It can be done by hand or automatically using Q and twoSz

%% Generating the goal representation based on the twoSz and Q
% It can be done by hand or automatically using Q and twoSz

if abs(twoSz)<eps
    REP_twoSz=1;
elseif twoSz > 0
    REP_twoSz = 2*twoSz + 1;
else
    REP_twoSz = - 2*twoSz;
end

if abs(Q)<eps
    REP_Q=1;
elseif Q > 0
    REP_Q = 2*Q + 1;
else
    REP_Q = - 2*Q;
end


%% Construct the goal representation based on the symmetries provided.
%% The order in the representation array must match the order of symmetries.

% Corresponding to symmetry_names ={U1_Q, U1_twoSz} in the input file.
if   no_of_symmetries == 2 && strcmp(symmetry_names{1}, 'U1_Q')  && strcmp(symmetry_names{2},'U1_twoSz')
    REP = [ REP_Q, REP_twoSz];
end





%% Generate the sites to be used

spinhalf_superfermion_site = ABSITE_spinhalf_superfermion(symmetry_names, convention);

%% Construct the  sites along the chain.

Sites = cell(1,L);
for i = 1:L
    Sites{i} = spinhalf_superfermion_site;
end



%% Generate the initial state. This is the infinite temperature state.
% The state corresponds to  |S> = |0, 0>  - i* |1,1>
%  at each site and corresponds to the trivial representation R=1;

cut = L;
stateMPS = ABMPS_create(L,cut,no_of_symmetries);
infMPS = ABMPS_create(L,cut,no_of_symmetries);
% if we use the convention='real' and U1_Q symmetry for the creation
% superoperators


%% construct the inf state. 

if   no_of_symmetries == 2 && strcmp(symmetry_names{1}, 'U1_Q') && strcmp(symmetry_names{2}, 'U1_twoSz')  && strcmp(convention, 'real')
    A = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},no_of_symmetries);    
    A = NTset_block(A,{{'t_in',1},{'tau',1},{'t_out',1}},{[1,1],[1,1],[1,1]},{'t_in','tau','t_out'},1/sqrt(4)*[1.0, 1.0, -1i, -1.i]);
   
end

if   no_of_symmetries == 2 && strcmp(symmetry_names{1}, 'U1_Q') && strcmp(symmetry_names{2}, 'U1_twoSz')  && strcmp(convention, 'real')
    schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},no_of_symmetries);
    schmidt = NTset_block(schmidt,{{'t_left',1},{'t_right',1}},{[1,1],[1,1],[1,1]},{'t_left','t_right'},[1.0]);
end



for site=1:L
   infMPS = ABMPS_set_matrix(infMPS,site,A); 
    infMPS.schmidt_list{site} = schmidt;
end

infMPS = ABMPS_set_schmidt(infMPS,schmidt);

%% Construct the initial infinite temperature state with some extra chemical potential.


if   no_of_symmetries == 2 && strcmp(symmetry_names{1}, 'U1_Q') && strcmp(symmetry_names{2}, 'U1_twoSz')  && strcmp(convention, 'real')
    B1 = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},no_of_symmetries);
    B2 = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},no_of_symmetries);
    
    vect = [1.0-mu ,1.0+mu, (-1.0)*1.0i , (-1.0)*1.0i];
    amp = norm(vect);
    B1 = NTset_block(B1,{{'t_in',1},{'tau',1},{'t_out',1}},{[1,1],[1,1],[1,1]},{'t_in','tau','t_out'},1/amp*vect);

    vect = [1.0+mu ,1.0-mu, (-1.0)*1.0i , (-1.0)*1.0i];
    amp = norm(vect);
    B2 = NTset_block(B2,{{'t_in',1},{'tau',1},{'t_out',1}},{[1,1],[1,1],[1,1]},{'t_in','tau','t_out'},1/amp*vect);

end


for site=1:L
    if site<=L/2
        stateMPS = ABMPS_set_matrix(stateMPS,site,B1);
    end
    
    if site >L/2
        stateMPS = ABMPS_set_matrix(stateMPS,site,B2);
    end
    
    stateMPS.schmidt_list{site} = schmidt;
end

schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},no_of_symmetries);
schmidt = NTset_block(schmidt,{{'t_left',1},{'t_right',1}},{[1,1],[1,1],[1,1]},{'t_left','t_right'},[1.0]);

stateMPS = ABMPS_set_schmidt(stateMPS,schmidt);

norm1= ABMPS_scalarprod(infMPS,stateMPS);
fprintf('Initial state norm = %.10f\n', real(norm1));


%% Initialize various variables of interest.
%occupation
average_occupation_TEBD = zeros(iteration_m, L);
%operators entanglement
Sop_TEBD = zeros(iteration_m, L-1);
%Renyi Entropy
Sre_TEBD = zeros(iteration_m);
%average current at the interface
average_current_TEBD = zeros(iteration_m,1);
average_current_u1_TEBD = zeros(iteration_m,1);
average_current_u2_TEBD = zeros(iteration_m,1);
average_current_d1_TEBD = zeros(iteration_m,1);
average_current_d2_TEBD = zeros(iteration_m,1);


%% Seting up the TEBD environment
fprintf ('Constructing the  TEBD environment...\n')

ABTEBD_env = ABTEBD_init(stateMPS,symmetries ,Sites,'FERMIONS');

%% Construct the Lindbladian and update the TEBD environment.
endbond = stateMPS.chain_length - 1;
fprintf ('Constructing the  Lindbladian...\n')

for bondpos = 1:endbond
    
    if  bondpos == 1
       
        ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env,bondpos,{...
            {'fdag_u', 'f_u', -0.5}, ...
            {'fdag_d', 'f_d', -0.5}, ...
            {'f_u', 'fdag_u', +0.5}, ...
            {'f_d', 'fdag_d', +0.5}, ...
            {'hub', 'id', 1.0*U}, ...
            {'id', 'hub', 0.5*U},...
            {'tilde_fdag_u', 'tilde_f_u', +0.5}, ...
            {'tilde_fdag_d', 'tilde_f_d', +0.5}, ...
            {'tilde_f_u', 'tilde_fdag_u', -0.5}, ...
            {'tilde_f_d', 'tilde_fdag_d', -0.5}, ...
            {'tilde_hub', 'id', -1.0*U},...
            {'id', 'tilde_hub', -0.5*U}...
            });
        
    elseif bondpos == L-1
        
        ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env,bondpos,{...
            {'fdag_u', 'f_u', -0.5}, ...
            {'fdag_d', 'f_d', -0.5}, ...
            {'f_u', 'fdag_u', +0.5}, ...
            {'f_d', 'fdag_d', +0.5}, ...
            {'id', 'hub', 1.0*U}, ...
            {'hub', 'id', 0.5*U}, ...
            {'tilde_fdag_u', 'tilde_f_u', +0.5}, ...
            {'tilde_fdag_d', 'tilde_f_d', +0.5}, ...
            {'tilde_f_u', 'tilde_fdag_u', -0.5}, ...
            {'tilde_f_d', 'tilde_fdag_d', -0.5}, ...
            {'tilde_hub', 'id', -0.5*U},...
            {'id', 'tilde_hub', -1.0*U}...
            });
    else
        ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env,bondpos,{...
            {'fdag_u', 'f_u', -0.5}, ...
            {'fdag_d', 'f_d', -0.5}, ...
            {'f_u', 'fdag_u', 0.5}, ...
            {'f_d', 'fdag_d', 0.5}, ...
            {'hub', 'id', 0.5*U}, ...
            {'id', 'hub', 0.5*U}, ...
            {'tilde_fdag_u', 'tilde_f_u', +0.5}, ...
            {'tilde_fdag_d', 'tilde_f_d', +0.5}, ...
            {'tilde_f_u', 'tilde_fdag_u', -0.5}, ...
            {'tilde_f_d', 'tilde_fdag_d', -0.5}, ...
            {'tilde_hub', 'id', -0.5*U}...
            {'id', 'tilde_hub', -0.5*U}...
            });
    end
end


%% Add the evolver to the TEBD environment.
fprintf ('Constructing the  evolver...\n')
ABTEBD_env = ABTEBD_gen_Ulist(ABTEBD_env, dt ,30);


%sites = [1:stateMPS.chain_length];
%occupation = average_occupation_TEBD(1,:);
%p = plot(sites,occupation, 'ro-');
%xlim([0 stateMPS.chain_length+1])
%ylim([0,1]);
%xlabel('sites');
%ylabel('<n>');

%p.XDataSource = 'sites';
%p.YDataSource = 'occupation';

%% Perform the TEBD evolution and evaluate the averages.

%index for the measurement step.
%performing the iterations
fprintf ('Performing the iterations...\n');

iteration_m = 0;
for iteration = 0:iteration_no-1
     tic
    if  mod(iteration, ratio)==0  || iteration == iteration_no-1
        
        iteration_m =iteration_m+1;
        time_m(iteration_m) = iteration*dt;
        
        % Calculate the norm at each iteration
        norm1= ABMPS_scalarprod(infMPS,ABTEBD_env.StateMPS);
        
        fprintf('State norm = %.10f\n', norm1);
        
        
        %compute the average occupation
        parfor site = 1:ABTEBD_env.chain_length
            stateMPS1 = ABMPS_apply_singlesite_op(ABTEBD_env.StateMPS,...
                site,...
                ABTEBD_env.Sites{site}.operators('n'),...
                [1,1],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                1);
            
            average_occupation_TEBD(iteration_m, site) = real(ABMPS_scalarprod(infMPS,stateMPS1)/ norm1);
        end
        
        
        fprintf ('N(%.2f) =',iteration*dt);
        display (average_occupation_TEBD(iteration_m, :));
        fprintf('N_tot(%.2f)= %.5f\n', iteration*dt, sum(average_occupation_TEBD(iteration_m, :)));
        occupation = average_occupation_TEBD(iteration_m, :);
       %refreshdata
       %drawnow
       %compute the average current across the interface
        
         site = ABTEBD_env.chain_length/2;
            
            %fdag_u*f_u
            stateMPS1 = ABMPS_apply_singlesite_op(ABTEBD_env.StateMPS,...
                site,...
                ABTEBD_env.Sites{site}.operators('fdag_u'),...
                [3,3],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            stateMPS1 = ABMPS_apply_singlesite_op(stateMPS1,...
                site+1,...
                ABTEBD_env.Sites{site+1}.operators('f_u'),...
                [2,2],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            %keyboard;
            average_current_u1_TEBD(iteration_m) =  -imag(ABMPS_scalarprod(infMPS,stateMPS1)/norm1);
            
            
            
            %-fu*fdag_u
            stateMPS1 = ABMPS_apply_singlesite_op(ABTEBD_env.StateMPS,...
                site,...
                ABTEBD_env.Sites{site}.operators('f_u'),...
                [2,2],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            stateMPS1 = ABMPS_apply_singlesite_op(stateMPS1,...
                site+1,...
                ABTEBD_env.Sites{site+1}.operators('fdag_u'),...
                [3,3],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            %keyboard;
            average_current_u2_TEBD(iteration_m) =  - imag(ABMPS_scalarprod(infMPS,stateMPS1)/norm1);
            
            
             %fdag_d*f_d
            stateMPS1 = ABMPS_apply_singlesite_op(ABTEBD_env.StateMPS,...
                site,...
                ABTEBD_env.Sites{site}.operators('fdag_d'),...
                [3,2],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            stateMPS1 = ABMPS_apply_singlesite_op(stateMPS1,...
                site+1,...
                ABTEBD_env.Sites{site+1}.operators('f_d'),...
                [2,3],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            %keyboard;
            average_current_d1_TEBD(iteration_m) =  -imag(ABMPS_scalarprod(infMPS,stateMPS1)/norm1);
            
            
            
             %f_d*fdag_d
            stateMPS1 = ABMPS_apply_singlesite_op(ABTEBD_env.StateMPS,...
                site,...
                ABTEBD_env.Sites{site}.operators('f_d'),...
                [2,3],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            stateMPS1 = ABMPS_apply_singlesite_op(stateMPS1,...
                site+1,...
                ABTEBD_env.Sites{site+1}.operators('fdag_d'),...
                [3,2],...
                ABTEBD_env.Sites,...
                ABTEBD_env.Symmetries,...
                -1);
            
            
            %keyboard;
            average_current_d2_TEBD(iteration_m) =  - imag(ABMPS_scalarprod(infMPS,stateMPS1)/norm1);
            
           
        
        fprintf ('J(%.2f) =',iteration*dt);
        
        average_current_TEBD(iteration_m)= average_current_u1_TEBD(iteration_m)+...
                                               average_current_u2_TEBD(iteration_m)+...  
                                               average_current_d1_TEBD(iteration_m)+...
                                               average_current_d2_TEBD(iteration_m);
                                           
        display(average_current_TEBD(iteration_m));
%         display(average_current_u1_TEBD(iteration_m, :));
%         display(average_current_u2_TEBD(iteration_m, :));
%         display(average_current_d1_TEBD(iteration_m, :));
%         display(average_current_d2_TEBD(iteration_m, :));
%         
        %compute the operators entanglement
        
%         [Sop, ~] = ABTEBD_get_entropy(ABTEBD_env);
%         
%         Sop_TEBD(iteration_m, :) = real(Sop);
%         
%         fprintf ('Sop(%.2f) =',iteration*dt);
%         display (Sop_TEBD(iteration_m, :));
        
        %         %Compute the Renyie entropy
        %
        %                 Sre_TEBD(iteration_m) = -log(ABMPS_scalarprod(ABTEBD_env.StateMPS,ABTEBD_env.StateMPS));
        %                 fprintf ('Sre(%.2f) =%.10f\n',iteration*dt, Sre_TEBD(iteration_m));
    end
    
    %perform one step in evolution
    [ABTEBD_env, info] = ABTEBD_evolve_Trotter1(ABTEBD_env, Mmax_TEBD, 1);
    fprintf ('M(%.2f) =',iteration*dt);
    display(info.bond_dims);
    fprintf ('Log(Truncerr(%.2f)) =',iteration*dt);
    display(real(log10(info.trunc_errs)));
    toc
    
end

save ('datafile.mat', 'average_occupation_TEBD', 'average_current_TEBD','Sop_TEBD');

return;

