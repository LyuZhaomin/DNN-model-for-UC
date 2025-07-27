% Step 0: Define parameters and plot a sample UC solution
define_constants;
cd(tempdir)                 % Go to a temporary writable folder
mpc = loadcase('case9');

N_gen = size(mpc.gen, 1);
N_bus = size(mpc.bus, 1);
% The ramp limit, start-up/shut-down cost, and 
% minimum up/down time are not specified in case118
mpc.gen(:, RAMP_10) = randi([90, 114], N_gen, 1);
mpc.gen(:, RAMP_30) = randi([90, 114], N_gen, 1);
mpc.gencost(:, STARTUP) = 1000 + (5000-1000)*rand(N_gen, 1);
mpc.gencost(:, SHUTDOWN) = 0.1*mpc.gencost(:, STARTUP);

savecase('case9_custom', mpc);  
mpc_test = loadcase('case9_custom');

baseMVA = mpc.baseMVA;
mui = randi([2,8], N_gen, 1);
mdi = max(1, mui - 2); 

mpopt = mpoption('model', 'DC', 'verbose', 0, 'out.all', 0);
results = runpf(mpc_test, mpopt);
PTDF = makePTDF(mpc_test);

T = 24;

URi = mpc_test.gen(:, RAMP_10);
DRi = mpc_test.gen(:, RAMP_30);
P_min = mpc_test.gen(:, PMIN);
P_max = mpc_test.gen(:, PMAX);

P_real = readmatrix('C:\Users\18216\Desktop\sample22.csv');
n_day = 5;
P_reshaped = reshape(P_real, 4, 24, n_day);  % [4 samples/hour × 24 hours × n_day]
P_hourly = squeeze(mean(P_reshaped, 1));  % [24 hours × n days]
P_mean_day = mean(P_hourly, 2);             
% 24×1 vector
load_shape = P_mean_day / mean(P_mean_day);  % normalized load shape (mean = 1)

base_demand = mpc_test.bus(:, PD);   
% N_bus × 1, the jth bus is valued if attached to demand.
                    
% Hinted by DeepOPF, the P_d here is randomly chosen from a sampling range.
% The range is [(1-x)base_demand, (1+x)base_demand], x could be 10%
% This is for input data generation. 
P_d = base_demand * load_shape'; % broadcast to N_bus × T
                                % load_shape' means transpose!

PF_ij_max = mpc_test.branch(:, RATE_A);
su_cost = mpc_test.gencost(:, 2); % (:,2) means the second column
sd_cost = mpc_test.gencost(:, 3);
a = mpc_test.gencost(:, COST);
b = mpc_test.gencost(:, COST+1);
c = mpc_test.gencost(:, COST+2);

G2B = zeros(N_bus, N_gen);
for i = 1:N_gen
    bus_index = mpc_test.gen(i, GEN_BUS);  % generator bus index
    G2B(bus_index, i) = 1;
end

PTDF_gen = PTDF * G2B;  % [N_branch × N_gen]
%% UC Formulation
% variables:
% sdpvar and binvar are used to define symbolic decision variables, with
% the syntax sdpvar(n, m, 'full') meaning a full n*m matrix.

SU = binvar(N_gen, T, 'full'); 
SD = binvar(N_gen, T, 'full');
O = binvar(N_gen, T, 'full');
P = sdpvar(N_gen, T, 'full');

% UC_obj = sum(sum(su_cost .* SU + sd_cost .* SD + a .* P.^2 + b .* P + c .* O));
UC_obj = sum(sum(repmat(su_cost,1,T) .* SU + repmat(sd_cost,1,T) .* SD + ...
                 repmat(a,1,T) .* P.^2 + repmat(b,1,T) .* P + repmat(c,1,T) .* O));

UC_constraints = {};  
for i = 1:N_gen
    for t = 2:T
        h_max_up = min(mui(i) + t - 1, T);
        for h = (t+1):h_max_up
            UC_constraints{end+1} = O(i,t) - O(i,t-1) - O(i,h) <= 0;
        end
        h_max_down = min(mdi(i) + t - 1, T);
        for h = (t+1):h_max_down
            UC_constraints{end+1} = O(i,t-1) - O(i,t) + O(i,h) <= 1;
        end
        UC_constraints{end+1} = O(i,t-1) - O(i,t) + SU(i,t) >= 0;
        UC_constraints{end+1} = O(i,t) - O(i,t-1) + SD(i,t) >= 0;
        UC_constraints{end+1} = P(i,t) - P(i,t-1) <= URi(i)*O(i,t-1) + ...
                                 P_min(i)*(O(i,t)-O(i,t-1)) + ...
                                 P_max(i)*(1 - O(i,t));
        UC_constraints{end+1} = P(i,t-1) - P(i,t) <= DRi(i)*O(i,t) + ...
                                 P_min(i)*(O(i,t-1)-O(i,t)) + ...
                                 P_max(i)*(1 - O(i,t-1));
    end
end

for i = 1:N_gen
    for t = 1:T
        UC_constraints{end+1} = P_min(i)*O(i,t) <= P(i,t) <= P_max(i)*O(i,t);
    end
end

for t = 1:T
    UC_constraints{end+1} = sum(P(:,t)) == sum(P_d(:,t));
end

for t = 1:T
    flow_t = PTDF_gen * P(:,t);  % branch flows due to generation
    for l = 1:size(mpc_test.branch,1)
        F_l = mpc_test.branch(l, RATE_A);  % thermal limit in MW
        if F_l == 0
            continue  % skip lines without limit
        end
        UC_constraints{end+1} = -F_l <= flow_t(l) <= F_l;
    end
end

UC_constraints = [UC_constraints{:}];
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
optimize(UC_constraints, UC_obj, options);

%% show the output, a combination of generator commitment schedule and generation dipatch.

O_value = value(O);   % N_gen × T
P_value = value(P);   % N_gen × T
%O_vec = reshape(O_value, 1, []);  % 1 × (N_gen*T)
%P_vec = reshape(P_value, 1, []);  % 1 × (N_gen*T)
Y_vec = [O_value, P_value];  % 1 × 144 (for N_gen=3, T=24)

%% visualization of the output data
figure;
plot(P_value', 'LineWidth', 1.5);
title('Generation Dispatch Over Time');
xlabel('Hour');
ylabel('Power Output (MW)');
legend({'Gen 1', 'Gen 2', 'Gen 3'});
grid on;

% generate the data for saving
% Step 1: Prepare for Data Generation
num_samples = 1000;  % Total data points
inputs = zeros(N_bus, T, num_samples);   % To store input (load)
outputs = zeros(N_gen, T, num_samples);  % To store output (generation)

% Step 2: Loop to Generate Data
sample_count = 0;
while sample_count < num_samples
    % 1. Generate new demand profile
    random_scale = 0.9 + 0.2 * rand(N_bus, 1);
    P_d = base_demand .* random_scale * load_shape';
    SU = binvar(N_gen, T, 'full'); 
    SD = binvar(N_gen, T, 'full');
    O = binvar(N_gen, T, 'full');
    P = sdpvar(N_gen, T, 'full');
    UC_obj = sum(sum(repmat(su_cost,1,T) .* SU + repmat(sd_cost,1,T) .* SD + ...
                 repmat(a,1,T) .* P.^2 + repmat(b,1,T) .* P + repmat(c,1,T) .* O));

    % 2. COMPLETELY reset constraints
    UC_constraints = {};
    
    % 3. REBUILD ALL CONSTRAINTS (in this order)
    
    % A. Power Balance (depends on P_d)
    for t = 1:T
        UC_constraints{end+1} = sum(P(:,t)) == sum(P_d(:,t));
    end
    
    % B. Line Flow (depends on P_d via PTDF)
    for t = 1:T
        flow_t = PTDF_gen * P(:,t);
        for l = 1:size(mpc_test.branch,1)
            if mpc_test.branch(l, RATE_A) > 0
                UC_constraints{end+1} = -mpc_test.branch(l, RATE_A) <= flow_t(l) <= mpc_test.branch(l, RATE_A);
            end
        end
    end
    
    % C. Generator Operational Constraints (MUST be rebuilt)
    for i = 1:N_gen
        for t = 2:T
            % Min up/down constraints
            h_max_up = min(mui(i) + t - 1, T);
            for h = (t+1):h_max_up
                UC_constraints{end+1} = O(i,t) - O(i,t-1) - O(i,h) <= 0;
            end
            h_max_down = min(mdi(i) + t - 1, T);
            for h = (t+1):h_max_down
                UC_constraints{end+1} = O(i,t-1) - O(i,t) + O(i,h) <= 1;
            end
            
            % Startup/shutdown logic
            UC_constraints{end+1} = O(i,t-1) - O(i,t) + SU(i,t) >= 0;
            UC_constraints{end+1} = O(i,t) - O(i,t-1) + SD(i,t) >= 0;
            
            % Ramp rate constraints
            UC_constraints{end+1} = P(i,t) - P(i,t-1) <= URi(i)*O(i,t-1) + ...
                                     P_min(i)*(O(i,t)-O(i,t-1)) + ...
                                     P_max(i)*(1 - O(i,t));
            UC_constraints{end+1} = P(i,t-1) - P(i,t) <= DRi(i)*O(i,t) + ...
                                     P_min(i)*(O(i,t-1)-O(i,t)) + ...
                                     P_max(i)*(1 - O(i,t-1));
        end
    end
    
    % D. Generation Limits
    for i = 1:N_gen
        for t = 1:T
            UC_constraints{end+1} = P_min(i)*O(i,t) <= P(i,t) <= P_max(i)*O(i,t);
        end
    end
    
    UC_constraints = [UC_constraints{:}];
    % 4. Solve fresh problem
    diagnostics = optimize(UC_constraints, UC_obj, options);
    
    % 5. Store if feasible
    if diagnostics.problem == 0
        sample_count = sample_count + 1;
        inputs(:,:,sample_count) = P_d;
        outputs(:,:,sample_count) = value(P);
    end
end
     
% Step 3: Save Dataset
save('UC_dataset_case9.mat', 'inputs', 'outputs');

% Step 4: Split into training/test sets
X_train = inputs(:, :, 1:800);
Y_train = outputs(:, :, 1:800);
X_test = inputs(:, :, 801:end);
Y_test = outputs(:, :, 801:end);

save('UC_train_test_case9.mat', 'X_train', 'Y_train', 'X_test', 'Y_test');
