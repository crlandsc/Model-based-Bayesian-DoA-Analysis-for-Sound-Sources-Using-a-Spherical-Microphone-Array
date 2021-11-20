%% BAYESIAN DoA ANALYSIS
clear all; close all; clc;
% Add Spherical Harmonics Functions
addpath('FUNCTIONS - SH and Beamforming');
% Add Experimental Data
addpath('Generated Noise Sources');

%% USER INPUTS
% =========================================================================
% ========================== INITIAL DEFINITIONS ==========================
% =========================================================================
% Enter location to save output folder. This will create a new folder in
% this location and save all figures and the workspace to this location
Active_Folder = 'C:\Users\Chris\Documents\Matlab\RPI Matlab\Research\Chris Research\';
% -------------------------------------------------------------------------
% Enter name of folder to save output files. This will create a new folder
% in this location and save all figures and the workspace to this location
Save_Folder_Name = '180705 - TEST';
Save_Location = strcat(Active_Folder,Save_Folder_Name);
mkdir(Save_Location)
% -------------------------------------------------------------------------
% Experimentally measured sources vs. artificially generated sources.
%       - Enter 1 to load experimentally measured sources
%       - Enter 0 to load artificially generated sources
Load_File_Yes_No = 0; % Enter 1 to load file, 0 to not load file
% -------------------------------------------------------------------------
% Enter name of file with experimentally measured data to be loaded
%       - (Files found in 'Generated Noise Sources' folder)
Load_Filename = 'WN_90_90_Windowed';
% -------------------------------------------------------------------------
% Nested sampling diffusion coefficient (Makes max less peaky and spreads
% out 'volume'). Recommended: Beta = 1.
Beta = 1; % 0 to 1 (0 = most diffuse, 1 = not diffuse)
% -------------------------------------------------------------------------
% Define size of initial table (exploration table).
%       - Recommended: 1000
Model_iterations = 1000;
% -------------------------------------------------------------------------
% Amount of times delta likelihood must be below Converge_Limit in order to
% break the loop (determines if delta converges)
Threshold_Amount = Model_iterations;
% -------------------------------------------------------------------------
% Define number of source models to run through.
%       - Recommended: 2 more than number of sources being tests. (e.g. If
%         3 sources are being tested, test up to models with 5 sources)
num_sources_guess = 5;
% -------------------------------------------------------------------------
% Define sound energy beam grid, a.k.a. resolution of ourput plots.
%       - Recommended: num_look_az = 100; num_look_el = 50;
num_look_az = 100; % Number of azimuthal look directions for grid
num_look_el = 50; % Number of elevation look directions for grid
% -------------------------------------------------------------------------
% Number of iterations before likelihood plot updates
%       *This does not affect any output, it just plots likelihood curve to
%       track progress*
Likelihood_Update = 10000;
% -------------------------------------------------------------------------
% Limit of times with no succesful step before alternative sampling is used
Stuck_Lim = 10; %0;
% -------------------------------------------------------------------------
% Number of times Likelihood increases before alternative sampling kicks in
L_up_idx_lim = 10000;
% -------------------------------------------------------------------------
% The threshold of delta each step must be below in order to "converge"
Converge_Limit = 10;
% -------------------------------------------------------------------------
% Threshold of the difference of max and min in the intial table in order
% to break
Exploration_Table_Threshold = 10;
% -------------------------------------------------------------------------


% =========================================================================
% ======================== CONSISTENT DEFINITIONS =========================
% =========================================================================
% User Input (Typically remains constant)
UI = struct; % Define User Input (UI) structure
% -------------------------------------------------------------------------
% Assign look direction grid for experimental model calcs
UI.num_look_az = num_look_az;
UI.num_look_el = num_look_el;
% -------------------------------------------------------------------------
% Choose which order SH beamformer is desired (order=2 for 16-channel mic)
order = 2;
UI.order = order;
% -------------------------------------------------------------------------
% Define boundary limits of look directions
phi_lim = [0 360];
UI.phi_lim = phi_lim;
theta_lim = [0 180];
UI.theta_lim = theta_lim;
% -------------------------------------------------------------------------
% Type of Beamformer
%       - Enter 1 for Plane Wave Decomposition Beamformer
%       - Enter 2 for Delay & Sum Beamformer
%       *All tests were run with plane wave decomp beamformer*
UI.Beam_Type = 1;
% -------------------------------------------------------------------------
% Calculate:
%       - Order index (ensures order (n) corresponds to appropriate degree (m))
%       - Look directions in spherical harmonics
%       - Modal amplitude in spherica harmonics
%       - Weights in spherical harmonics
[UI.order_idx,UI.ldir_harms,UI.b_n,UI.W_nm] = Beam_Model_Precalcs_4_7_18(UI);
UI.plotflag = 0; % supress plots for now (don't change)

% =========================================================================
% =============== ARTIFICIALLY GENERATED SOURCE ATTRIBUTES ================
%     ***Ignore this section if testing experimentally measured data***
% =========================================================================
if Load_File_Yes_No == 0
% Artificially generated DoA attributes. Define data to be tested against.
UI.amplitudes = [1 1 1]; % 0 to 1 amplitude(s) per source(s)
UI.S_angles_az = [5 135 270]; % 0 to 360 (degrees) azimuthal location(s) for source(s)
UI.S_angles_el = [60 140 90]; % 0 to 180 (degrees) elevation location(s) for source(s)
% -------------------------------------------------------------------------
% Amplitude of background noise.
%       - Recommended: between 0 and .3
min_Noise_Level = .05;
max_Noise_Level = .25;
Noise_Level = max_Noise_Level - min_Noise_Level;
% -------------------------------------------------------------------------
% Number of periods for noise (Changes behavior of noise)
x_periods = 1;
y_periods = 1;
% -------------------------------------------------------------------------
% Generate Noise
[Noise] = Beam_Noise(UI,min_Noise_Level,max_Noise_Level);
% -------------------------------------------------------------------------
% Generate Artificial Sources
[Beam_Exp_No_Noise] = Beam_Model_Iterate_4_7_18(UI); % Beamform artificial data
% -------------------------------------------------------------------------
% Add noise to artificially generated signals.
% Weight noise to affect low signals more & beam locations less.
Inv_Beam_Exp_No_Noise = 1-Beam_Exp_No_Noise; % Inverse generated beam grid
Beam_Exp = Beam_Exp_No_Noise + (Inv_Beam_Exp_No_Noise.^5.*Noise); % Add noise
Beam_Exp = Beam_Exp/ max(max(Beam_Exp)); % Normalize
% -------------------------------------------------------------------------
% Plot & save noise and artificially generated beam source data
Plot_Beam_Grid(Noise,UI);
savefig(strcat(Save_Location,'\Noise.fig'));
Plot_Beam_Grid(Beam_Exp,UI);
savefig(strcat(Save_Location,'\Beam_Exp.fig'));
% -------------------------------------------------------------------------
% Calculate evidence of single source model against background noise (e.g. no sources)
Background_Noise_Evidence = Beam_Noise_Evidence_New(UI,Noise,Model_iterations);


% =========================================================================
% ============== EXPERIMENTALLY MEASURED SOURCE ATTRIBUTES ================
else %***Ignore this section if testing artificially generated data***
% =========================================================================
% Load experimental recorded data
load(strcat(Load_Filename,'.mat')); % Input signal
UI.fs = fs; % Sampling frequency
UI.mic_signals = WN; % Data is saved as 'WN' in each file.
% -------------------------------------------------------------------------
% Calculate, plot, & save experimentally measured beam source data
UI.plotflag = 1; % plot experimental data in order to save
[Beam_Exp] = Beam_Data_4_7_18(UI); % Beamform experimental data
savefig(strcat(Save_Location,'\Beam_Exp.fig'));
% -------------------------------------------------------------------------
% Evidence of single source model against background noise (e.g. no sources)
% Table = 1000; phi x theta = 100x50; beta = 1; 
Background_Noise_Evidence = -6086.65229863351;
end


%% INITIAL DEFINITIONS
% ============================ Initializations ============================
% These values will only be defined once before all iterations are run.
Amps = zeros(num_sources_guess,num_sources_guess);
Phis = Amps;
Thetas = Amps;
Sort_Amp_Std = Amps;
Sort_Phi_Std = Amps;
Sort_Theta_Std = Amps;

Completely_Sorted_Table = zeros(Model_iterations,16,num_sources_guess);

Shift_Final_L = zeros(1,num_sources_guess);
Parameter_Break_VS_Likelihood_Break = Shift_Final_L;
Z = Shift_Final_L;
Evidence_log_noshift = Shift_Final_L;
Evidence_log = Shift_Final_L;

count_stuck = zeros(1,num_sources_guess);
Direct_Slice_Sample_idx = zeros(1,num_sources_guess);
Random_Gen_Slice_Samp_idx = zeros(1,num_sources_guess);
B_i_j = zeros(1,num_sources_guess);
deciban = zeros(1,num_sources_guess);

Total_Time = 0;
idx_total = 0;

Converge_Limit_all = Converge_Limit;
Initial_Table_Threshold_all = Exploration_Table_Threshold;


%% MODEL SELECTION
% =========================================================================
% ---------------------- Model Selection Iterations -----------------------
% =========================================================================
% Begin iterations through model selections. Begin with a single sources
% model, then double source, and so on and so forth.
for source_idx = 1:num_sources_guess

if source_idx == 1
    Converge_Limit = 1;
    Exploration_Table_Threshold = 1;
else
    Converge_Limit = Converge_Limit_all;
    Exploration_Table_Threshold = Initial_Table_Threshold_all;
end

%% --------------------------- Initializations ----------------------------

% Calculate randomized DoA attributes to populate exploration table
theta = 180*rand(Model_iterations,source_idx);
phi = 360*rand(Model_iterations,source_idx);
amplitudes = [ones(Model_iterations,1) rand(Model_iterations,source_idx-1)];

% Initialize for Explorations Table and initial iterative calculations
iteration = 0; % Initialize iteration count
Beam_Model = zeros(UI.num_look_el,UI.num_look_az); % Initialize beam model grid
error = zeros(Model_iterations,1); % Initialize error
E = zeros(Model_iterations,1); % Initialize evidence
Likelihood = zeros(Model_iterations,1); % Initialize likelihood
K = UI.num_look_az * UI.num_look_el; % Calculate K (used in likelihood calculations, aka gridsize))
Exploration_Table = zeros(Model_iterations, 1+3*source_idx); % Initialize exploration table

% Stop Matlab from breaking (ensure there aren't a bajillion plots)
UI.plotflag = 0; % ***NEVER COMMENT OUT***

%% Calculate Initial Likelihoods and Populate Exploration Table
for idx = 1:Model_iterations
    
    % DoA Attributes
    UI.amplitudes = amplitudes(idx,:); % assign amplitude per iteration
    UI.S_angles_az = phi(idx,:); % assign source phi per iteration
    UI.S_angles_el = theta(idx,:); % assign source theta per iteration
    
    % Beamform (create model)
    Beam_Model = Beam_Model_Iterate_4_7_18(UI);
    
    % Calculate Likelihood from model
    difference = Beam_Exp - Beam_Model;
    error = sum(sum(difference.*difference));
    E = error/2;
    Likelihood(idx) = (-Beta*K/2)*log10(E);
    
    % Populate Exploration Table
    % Format: [Likelihood, Amplitudes, Phis, Thetas]
    Exploration_Table(idx,:) = [Likelihood(idx), UI.amplitudes, UI.S_angles_az, UI.S_angles_el];
end

% ------------------------- Model Initializations -------------------------
% These values will be reset to their original definitions for each
% different model iteration. They will not be reset for each parameter
% iteration
New_Likelihood = Likelihood;
New_amplitudes = amplitudes;
New_phi = phi;
New_theta = theta;

Sort_Amp = zeros(Model_iterations,source_idx);
Sort_Phi = zeros(Model_iterations,source_idx);
Sort_Theta = zeros(Model_iterations,source_idx);
std_max = zeros(1,3);
std_min = zeros(1,3);
New_Table = zeros(1,size(Exploration_Table,2));
Loop_Break = ones(1,Threshold_Amount);

Iteration_Time = 0;
Iteration_Time_1000 = 0;
idx = 0;
L_up_idx = 0;
Paramter_Dif_Break = 0;
stuck = 0;
check_stuck = 0;
Random_Gen_Slice_Samp = 0;

idk = 1;
std_idx = 1;
Run_plot_idx = 1;
Area_Converge = 1;
Use_Other_Sampling_Methods = 1;
Int = 1;

Delta_Like = inf;
Loop_Sum = inf;
Sort_Amp_Std_Avg = inf;
Sort_Phi_Std_Avg = inf;
Sort_Theta_Std_Avg = inf;
Check_Threshold = inf;
Max_L = max(Exploration_Table(:,1));


%% Loop until likelihood increase plateaus and converges
while Area_Converge == 1 && Paramter_Dif_Break ~= 1 || Loop_Sum >= 1 && Paramter_Dif_Break ~= 1
    tic;
    
    % Maximum likelihood should never eclipse zero
    if Max_L > 0
        break
    end
    
    % Count the amount of times no successful step is made (for slice sampling)
    stuck = stuck + 1;
    
    % Count total iterations
    idx = idx + 1;
    
    % Find max/min likehoods
    [Max_L,~] = max(New_Likelihood); % Minimum likelihood
    [Min_L,Min_idx_L] = min(New_Likelihood); % Minimum likelihood
    
% -------------------------- Random Perturbation --------------------------
    % Select lowest likelihood in the initial table of paramters and
    % randomly perturb values until the likelihood improves. If it
    % improves, replace values in the initial table and append them in a
    % new, ever expnding, table as well.
    if L_up_idx < L_up_idx_lim && Use_Other_Sampling_Methods == 1 && stuck < Stuck_Lim
        % Randomly select parameter(s)
        if source_idx == 1
            random = round((7.49999-1.5)*rand(1,1)+1.5);
        else
            random = round((7.49999-.5)*rand(1,1)+.5);
        end
        
        % Perturb random parameter(s)
        switch random
            case 1 % Perturb Amplitude
                slot_a = round((source_idx-1)*rand(1,1)+1);
                slot = round(source_idx*rand(1,1));
                if slot_a == 1
                    slot_a = source_idx;
                end
                Int_Amp = New_amplitudes(Min_idx_L,:); % assign amplitude to rerun calcs
                Int_Amp(:,slot_a) = rand(1,1); % insert perturbation
                Int_Phi = New_phi(Min_idx_L,:); % assign source phi to rerun calcs
                Int_Theta = New_theta(Min_idx_L,:); % assign source theta to rerun calcs
            case 2 % Perturb Phi
                slot_a = round((source_idx-1)*rand(1,1)+1);
                slot = round(source_idx*rand(1,1));
                if slot == 0
                    slot = source_idx;
                end
                Int_Amp = New_amplitudes(Min_idx_L,:); % assign amplitude to rerun calcs
                Int_Phi = New_phi(Min_idx_L,:); % assign source phi to rerun calcs
                Int_Phi(:,slot) = 360*rand(1,1); % insert perturbation
                Int_Theta = New_theta(Min_idx_L,:); % assign source theta to rerun calcs
            case 3 % Perturb Theta
                slot_a = round((source_idx-1)*rand(1,1)+1);
                slot = round(source_idx*rand(1,1));
                if slot == 0
                    slot = source_idx;
                end
                Int_Amp = New_amplitudes(Min_idx_L,:); % assign amplitude to rerun calcs
                Int_Phi = New_phi(Min_idx_L,:); % assign source phi to rerun calcs
                Int_Theta = New_theta(Min_idx_L,:); % assign source theta to rerun calcs
                Int_Theta(:,slot) = 180*rand(1,1); % insert perturbation
            case 4 % Perturb Amplitude & Phi
                slot_a = round((source_idx-1)*rand(1,1)+1);
                slot = round(source_idx*rand(1,1));
                if slot_a == 1
                    slot_a = source_idx;
                end
                if slot == 0
                    slot = source_idx;
                end
                Int_Amp = New_amplitudes(Min_idx_L,:); % assign amplitude to rerun calcs
                Int_Amp(:,slot_a) = rand(1,1); % insert perturbation
                Int_Phi = New_phi(Min_idx_L,:); % assign source phi to rerun calcs
                Int_Phi(:,slot) = 360*rand(1,1); % insert perturbation
                Int_Theta = New_theta(Min_idx_L,:); % assign source theta to rerun calcs
            case 5 % Perturb Amplitude & Theta
                slot_a = round((source_idx-1)*rand(1,1)+1);
                slot = round(source_idx*rand(1,1));
                if slot_a == 1
                    slot_a = source_idx;
                end
                if slot == 0
                    slot = source_idx;
                end
                Int_Amp = New_amplitudes(Min_idx_L,:); % assign amplitude to rerun calcs
                Int_Amp(:,slot_a) = rand(1,1); % insert perturbation
                Int_Phi = New_phi(Min_idx_L,:); % assign source phi to rerun calcs
                Int_Theta = New_theta(Min_idx_L,:); % assign source theta to rerun calcs
                Int_Theta(:,slot) = 180*rand(1,1); % insert perturbation
            case 6 % Perturb Phi & Theta
                slot_a = round((source_idx-1)*rand(1,1)+1);
                slot = round(source_idx*rand(1,1));
                if slot_a == 1
                    slot_a = source_idx;
                end
                if slot == 0
                    slot = source_idx;
                end
                Int_Amp = New_amplitudes(Min_idx_L,:); % assign amplitude to rerun calcs
                Int_Phi = New_phi(Min_idx_L,:); % assign source phi to rerun calcs
                Int_Phi(:,slot) = 360*rand(1,1); % insert perturbation
                Int_Theta = New_theta(Min_idx_L,:); % assign source theta to rerun calcs
                Int_Theta(:,slot) = 180*rand(1,1); % insert perturbation
            case 7 % Perturb Amplitude, Phi, & Theta
                slot_a = round((source_idx-1)*rand(1,1)+1);
                slot = round(source_idx*rand(1,1));
                if slot_a == 1
                    slot_a = source_idx;
                end
                if slot == 0
                    slot = source_idx;
                end
                Int_Amp = New_amplitudes(Min_idx_L,:); % assign amplitude to rerun calcs
                Int_Amp(:,slot_a) = rand(1,1); % insert perturbation
                Int_Phi = New_phi(Min_idx_L,:); % assign source phi to rerun calcs
                Int_Phi(:,slot) = 360*rand(1,1); % insert perturbation
                Int_Theta = New_theta(Min_idx_L,:); % assign source theta to rerun calcs
                Int_Theta(:,slot) = 180*rand(1,1); % insert perturbation
        end % end switch random
        check_stuck = 0; % set stuck value count to 0 (don't need to go to slice sampling yet)
        
        
% -------------------- Alternative Random Perturbation --------------------
    elseif stuck >= Stuck_Lim && Use_Other_Sampling_Methods == 1 ||...
            Use_Other_Sampling_Methods == 1 && L_up_idx >= L_up_idx_lim
        
        % Randomly select parameter(s)
        if source_idx == 1
            random = round((7.49999-1.5)*rand(1,1)+1.5);
        else
            random = round((7.49999-.5)*rand(1,1)+.5);
        end
        
        % Sort exploration table by likelihood (least to greatest)
        Sort = sortrows(Exploration_Table,1);
        Sort_Parameters = Sort(:,2:end);
        
        % Select random model in exploration table for perturbation limit
        random_sort_idx = round((Model_iterations - 1)*rand(1) + 1);
        
        % Calculate the delta limits for perturbation limit per parameter type
        Sort_Delta_Amp = abs(sort(Sort_Parameters(1,1:source_idx)) -...
            sort(Sort_Parameters(random_sort_idx,1:source_idx)));
        Sort_Delta_Phi = abs(sort(Sort_Parameters(1,source_idx+1:2*source_idx)) -...
            sort(Sort_Parameters(random_sort_idx,source_idx+1:2*source_idx)));
        Sort_Delta_Theta = abs(sort(Sort_Parameters(1,2*source_idx+1:3*source_idx)) -...
            sort(Sort_Parameters(random_sort_idx,2*source_idx+1:3*source_idx)));
        
        % Calculate maximum per parameter type
        Max_Delta_Amp = max(Sort_Delta_Amp);
        Max_Delta_Phi = max(Sort_Delta_Phi);
        Max_Delta_Theta = max(Sort_Delta_Theta);
        
        % Perturb random parameter(s)
        switch random
            case 1 % Perturb Amplitude
                Delta_Amp = (2*Max_Delta_Amp)*rand(1,source_idx-1) - Max_Delta_Amp; % Delta Amp
                Int_Amp = [1 Sort_Parameters(1,2:source_idx)+Delta_Amp]; % Adjust value(s) by delta
                for Int_idx = 1:source_idx % Account for value(s) out of bounds
                    if Int_Amp(Int_idx) > 1
                        Int_Amp(Int_idx) = 1 - (Int_Amp(Int_idx)-1);
                    end
                end
                Int_Phi = Sort_Parameters(1,source_idx+1:2*source_idx); % Assign unaffected value(s)
                Int_Theta = Sort_Parameters(1,2*source_idx+1:3*source_idx); % Assign unaffected value(s)
            case 2 % Perturb Phi
                Delta_Phi = (2*Max_Delta_Phi)*rand(1,source_idx) - Max_Delta_Phi; % Delta Phi
                Int_Phi = Sort_Parameters(1,source_idx+1:2*source_idx) + Delta_Phi; % Adjust value(s) by delta
                for Int_idx = 1:source_idx % Account for value(s) out of bounds
                    if Int_Phi(Int_idx) > 360
                        Int_Phi(Int_idx) = Int_Phi(Int_idx) - 360;
                    end
                end
                Int_Amp = [1 Sort_Parameters(1,2:source_idx)]; % Assign unaffected value(s)
                Int_Theta = Sort_Parameters(1,2*source_idx+1:3*source_idx); % Assign unaffected value(s)
            case 3 % Perturb Theta
                Delta_Theta = (2*Max_Delta_Theta)*rand(1,source_idx) - Max_Delta_Theta; % Delta Theta
                Int_Theta = Sort_Parameters(1,2*source_idx+1:3*source_idx) + Delta_Theta; % Adjust value(s) by delta
                for Int_idx = 1:source_idx % Account for value(s) out of bounds
                    if Int_Theta(Int_idx) > 180
                        Int_Theta(Int_idx) = Int_Theta(Int_idx) - 360;
                    end
                end
                Int_Amp = [1 Sort_Parameters(1,2:source_idx)]; % Assign unaffected value(s)
                Int_Phi = Sort_Parameters(1,source_idx+1:2*source_idx); % Assign unaffected value(s)
            case 4 % Perturb Amplitude & Phi
                Delta_Amp = (2*Max_Delta_Amp)*rand(1,source_idx-1) - Max_Delta_Amp; % Delta Amp
                Int_Amp = [1 Sort_Parameters(1,2:source_idx) + Delta_Amp]; % Adjust value(s) by delta
                Delta_Phi = (2*Max_Delta_Phi)*rand(1,source_idx) - Max_Delta_Phi; % Delta Phi
                Int_Phi = Sort_Parameters(1,source_idx+1:2*source_idx) + Delta_Phi; % Adjust value(s) by delta
                for Int_idx = 1:source_idx % Account for value(s) out of bounds
                    if Int_Amp(Int_idx) > 1
                        Int_Amp(Int_idx) = 1 - (Int_Amp(Int_idx)-1);
                    end
                    if Int_Phi(Int_idx) > 360
                        Int_Phi(Int_idx) = Int_Phi(Int_idx) - 360;
                    end
                end
                Int_Theta = Sort_Parameters(1,2*source_idx+1:3*source_idx); % Assign unaffected value(s)
            case 5 % Perturb Amplitude & Theta
                Delta_Amp = (2*Max_Delta_Amp)*rand(1,source_idx-1) - Max_Delta_Amp; % Delta Amp
                Int_Amp = [1 Sort_Parameters(1,2:source_idx) + Delta_Amp]; % Adjust value(s) by delta
                Delta_Theta = (2*Max_Delta_Theta)*rand(1,source_idx) - Max_Delta_Theta; % Delta Theta
                Int_Theta = Sort_Parameters(1,2*source_idx+1:3*source_idx) + Delta_Theta; % Adjust value(s) by delta
                for Int_idx = 1:source_idx % Account for value(s) out of bounds
                    if Int_Amp(Int_idx) > 1
                        Int_Amp(Int_idx) = 1 - (Int_Amp(Int_idx)-1);
                    end
                    if Int_Theta(Int_idx) > 180
                        Int_Theta(Int_idx) = Int_Theta(Int_idx) - 180;
                    end
                end
                Int_Phi = Sort_Parameters(1,source_idx+1:2*source_idx); % Assign unaffected value(s)
            case 6 % Perturb Phi & Theta
                Delta_Phi = (2*Max_Delta_Phi)*rand(1,source_idx) - Max_Delta_Phi; % Delta Phi
                Int_Phi = Sort_Parameters(1,source_idx+1:2*source_idx) + Delta_Phi; % Adjust value(s) by delta
                Delta_Theta = (2*Max_Delta_Theta)*rand(1,source_idx) - Max_Delta_Theta; % Delta Theta
                Int_Theta = Sort_Parameters(1,2*source_idx+1:3*source_idx) + Delta_Theta; % Adjust value(s) by delta
                for Int_idx = 1:source_idx % Account for value(s) out of bounds
                    if Int_Phi(Int_idx) > 360
                        Int_Phi(Int_idx) = Int_Phi(Int_idx) - 360;
                    end
                    if Int_Theta(Int_idx) > 180
                        Int_Theta(Int_idx) = Int_Theta(Int_idx) -  180;
                    end
                end
                Int_Amp = [1 Sort_Parameters(1,2:source_idx)]; % Assign unaffected value(s)
            case 7 % Perturb Amplitude, Phi, & Theta
                Delta_Amp = (2*Max_Delta_Amp)*rand(1,source_idx-1) - Max_Delta_Amp; % Delta Amp
                Int_Amp = [1 Sort_Parameters(1,2:source_idx) + Delta_Amp]; % Adjust value(s) by delta
                Delta_Phi = (2*Max_Delta_Phi)*rand(1,source_idx) - Max_Delta_Phi; % Delta Phi
                Int_Phi = Sort_Parameters(1,source_idx+1:2*source_idx) + Delta_Phi; % Adjust value(s) by delta
                Delta_Theta = (2*Max_Delta_Theta)*rand(1,source_idx) - Max_Delta_Theta; % Delta Theta
                Int_Theta = Sort_Parameters(1,2*source_idx+1:3*source_idx) + Delta_Theta; % Adjust value(s) by delta
                for Int_idx = 1:source_idx % Account for value(s) out of bounds
                    if Int_Amp(Int_idx) > 1
                        Int_Amp(Int_idx) = 1 - (Int_Amp(Int_idx)-1);
                    end
                    if Int_Phi(Int_idx) > 360
                        Int_Phi(Int_idx) = Int_Phi(Int_idx) - 360;
                    end
                    if Int_Theta(Int_idx) > 180
                        Int_Theta(Int_idx) = Int_Theta(Int_idx) - 180;
                    end
                end
        end % end switch random
        
        % Prevent loop from getting stuck on one value.
        % Direct slice sampling (This should rarely happen in this context)
        if stuck > 20000
            rand_stuck = round((Model_iterations) * rand(1,1) + .5);
            if rand_stuck > Model_iterations
                rand_stuck = Model_iterations;
            end
            Int_Amp = Sort(rand_stuck,2:source_idx+1);
            Int_Phi = Sort(rand_stuck,source_idx+2:2*source_idx+1);
            Int_Theta = Sort(rand_stuck,2*source_idx+2:3*source_idx+1);
            Completely_Stuck_Check = 1;
        end
        
        check_stuck = 1; % notification of alternative sampling
    end
    
    % Assign Sampled Values
    UI.amplitudes = abs(Int_Amp);
    UI.S_angles_az = abs(Int_Phi);
    UI.S_angles_el = abs(Int_Theta);
    
% -------------------- Beamforming & Likelihood Calcs ---------------------
    % Beamform (create model)
    Beam_Model = Beam_Model_Iterate_4_7_18(UI);
    
    % Calculate Likelihood
    difference = Beam_Exp - Beam_Model;
    error = sum(sum(difference.*difference));
    E = error/2;
    Check_Likelihood = (-Beta*K/2)*log10(E);
    
    % Replace explorations table values
    % [Likelihood, Amplitudes, Phis, Thetas]
    % Dimensions: (Model_iterations) x (3*source_idx + 1)
    if Check_Likelihood > Min_L

        New_amplitudes(Min_idx_L,:) = UI.amplitudes; % assign new amplitude
        New_phi(Min_idx_L,:) = UI.S_angles_az; % assign new amplitude
        New_theta(Min_idx_L,:) = UI.S_angles_el; % assign new amplitude
        
        % Populate and Rewrite Parameter Tables
        L_up_idx = L_up_idx+1;
        New_Table(L_up_idx,:) = Exploration_Table(Min_idx_L,:); % Append unperturbed entry to to new table
        New_Likelihood(Min_idx_L) = Check_Likelihood;
        % Replace old entry in exploration table with higher likelihood entry
        Exploration_Table(Min_idx_L,:) =...
            [New_Likelihood(Min_idx_L), UI.amplitudes, UI.S_angles_az, UI.S_angles_el];
        
        % Calculate the delta of likelihood
        if L_up_idx > 1
            Delta_Like(L_up_idx) = New_Table(L_up_idx,1) - New_Table(L_up_idx-1,1);
        end
        
% --------------------------- Threshold Check -----------------------------
        % Check to see if likelihood has plateaued. 1 for no, 0 for yes.
        if L_up_idx > Threshold_Amount % Needs to have minimum # of entries
            Check_Threshold = abs(sum(Delta_Like(L_up_idx-Threshold_Amount+1:L_up_idx)));
            if Check_Threshold > Converge_Limit % If > convergence criteria, add 1
                Loop_Break = circshift(Loop_Break,-1);
                Loop_Break(Threshold_Amount) = 1;
            else                                % If < convergence criteria, add 0
                Loop_Break = circshift(Loop_Break,-1);
                Loop_Break(Threshold_Amount) = 0;
            end
            
            % Keep track of amount of times various sampling methods have
            % been used
            if check_stuck == 1
                count_stuck(source_idx) = count_stuck(source_idx) + 1;
            end
            if Completely_Stuck_Check == 1
                Direct_Slice_Sample_idx(source_idx) = Direct_Slice_Sample_idx(source_idx) + 1;
            end
            if Random_Gen_Slice_Samp == 1
                Random_Gen_Slice_Samp_idx(source_idx) = Random_Gen_Slice_Samp_idx(source_idx) + 1;
            end
            
            Loop_Sum = sum(Loop_Break); % Add together to determine if convergence criteria is met
        end
        stuck = 0; % Reset amount of times no successful step was made
    end
    
% ------------------------- Parameter Convergence -------------------------
    % Check differences between Paramters within table (ensure convergence)
    if idx == 5000*std_idx && source_idx > 1 % only every 5000 iterations
        for Table_idx = 1:Model_iterations            
            Sort_Amp(Table_idx,:) = sort(Exploration_Table(Table_idx,2:source_idx+1));
            Sort_Phi(Table_idx,:) = sort(Exploration_Table(Table_idx,source_idx+2:2*source_idx+1));
            Sort_Theta(Table_idx,:) = sort(Exploration_Table(Table_idx,2*source_idx+2:3*source_idx+1));
            Sort_Table_Low_to_High = [New_Likelihood Sort_Amp Sort_Phi Sort_Theta];
            
             % Organize Table for each model iteration for convenient
             % viewing if desired (not necessary for functionality of code)
            Completely_Sorted_Table(:,1:(1+(3*source_idx)),source_idx) = sortrows(Sort_Table_Low_to_High,1);
        end
        
        % Calculate 2*sigma to determine spread (95%) of parameters
        for Sort_Std_idx = 1:source_idx
            Sort_Amp_Std(source_idx,Sort_Std_idx) = 2.*std(Sort_Amp(:,Sort_Std_idx));
            Sort_Phi_Std(source_idx,Sort_Std_idx) = 2.*std(Sort_Phi(:,Sort_Std_idx));
            Sort_Theta_Std(source_idx,Sort_Std_idx) = 2.*std(Sort_Theta(:,Sort_Std_idx));
        end
        
        % Average the 2*sigma (in case there are multiple sources)
        Sort_Amp_Std_Avg = sum(Sort_Amp_Std(source_idx,:))/source_idx;
        Sort_Phi_Std_Avg = sum(Sort_Phi_Std(source_idx,:))/source_idx;
        Sort_Theta_Std_Avg = sum(Sort_Theta_Std(source_idx,:))/source_idx;
        
        % If parameters meet the convergence criteria, then break
        if Sort_Amp_Std_Avg < .01 && Sort_Phi_Std_Avg < 1 && Sort_Theta_Std_Avg < .5
            Paramter_Dif_Break = 1;
            % If this is 0, then the likelihood converged.
            % If this is 1, then the parameter became very similar, so
            % there was no room for improvement, thus converging
            Parameter_Break_VS_Likelihood_Break(source_idx) = 1;
            break
        end
        
        std_idx = std_idx + 1; % count number of times parameter convergence checked
    end % End if check for parameter convergence
    
% --------- Determine if there is still more area to be explored ----------
    % Calculate the range of likelihood in the exploration table
    L_Dif_Initial_Table = max(Exploration_Table(:,1)) - min(Exploration_Table(:,1));
    if L_Dif_Initial_Table >= Exploration_Table_Threshold
        Area_Converge = 1;
    else
        Area_Converge = 0;
    end
    
% --------------------------- Slice Sampling 1 ----------------------------
    % Completely repopulate exploration table & see if likelihood increases
    if Area_Converge == 1 && stuck == 200 && Direct_Slice_Sample_idx(source_idx) < 200
        
        Slice_amplitudes = [ones(Model_iterations,1) rand(Model_iterations,source_idx-1)];
        Slice_phi = 360*rand(Model_iterations,source_idx);
        Slice_theta = 180*rand(Model_iterations,source_idx);

        % Initialize for exploration Table and Iterative Calculations
        Slice_Table = zeros(Model_iterations, 1+3*source_idx);

        for Slice_idx = 1:Model_iterations

            % Iterative Calculations
            iteration = iteration + 1; % count iteration
            Int_Amp = Slice_amplitudes(Slice_idx,:); % assign amplitude per iteration
            UI.amplitudes = Int_Amp;
            Int_Phi = Slice_phi(Slice_idx,:); % assign source phi per iteration
            UI.S_angles_az = Int_Phi;
            Int_Theta = Slice_theta(Slice_idx,:); % assign source theta per iteration
            UI.S_angles_el = Int_Theta;

            % Beamform (create model)
            Beam_Model = Beam_Model_Iterate_4_7_18(UI);

            % Calculate Likelihood
            difference = Beam_Exp - Beam_Model;
            error = sum(sum(difference.*difference));
            E = error/2;
            Slice_Like = (-Beta*K/2)*log10(E);

            % [Likelihood, Amplitudes, Phis, Thetas]
            Slice_Table(Slice_idx,:) = [Slice_Like, Int_Amp, Int_Phi, Int_Theta];
        end

        [~,Max_Slice_Like_idx] = max(Slice_Table(:,1));
        Int_Amp = Slice_Table(Max_Slice_Like_idx,2:source_idx+1);
        Int_Phi = Slice_Table(Max_Slice_Like_idx,source_idx+2:2*source_idx+1);
        Int_Theta = Slice_Table(Max_Slice_Like_idx,2*source_idx+2:3*source_idx+1);
        Random_Gen_Slice_Samp = 1;
        Completely_Stuck_Check = 0;
        Use_Other_Sampling_Methods = 0;
        
% --------------------------- Slice Sampling 2 ----------------------------
    % Direct slice sampling
    elseif Area_Converge == 1 && stuck > 200 || Loop_Sum == 0 && Area_Converge == 1
        % Prevent loop from getting stuck on one value
        rand_stuck = round((Model_iterations) * rand(1,1) + .5);
        Sort = sortrows(Exploration_Table,1); % Sort initial table by likelihood (least to greatest)
        Int_Amp = Sort(rand_stuck,2:source_idx+1);
        Int_Phi = Sort(rand_stuck,source_idx+2:2*source_idx+1);
        Int_Theta = Sort(rand_stuck,2*source_idx+2:3*source_idx+1);
        Completely_Stuck_Check = 1;
        Random_Gen_Slice_Samp = 0;
        Use_Other_Sampling_Methods = 0;
    else
        %Otherwise use normal sampling methods
        Completely_Stuck_Check = 0;
        Random_Gen_Slice_Samp = 0;
        Use_Other_Sampling_Methods = 1;
    end
% ----------------------- End Sampling Calculations -----------------------
    
% ---------------------------- Progress Display ---------------------------
    % Define time & iterations
    single_time = toc;
    Iteration_Time = single_time + Iteration_Time;
    Iteration_Time_1000 = single_time + Iteration_Time_1000;
    Total_Time = Total_Time + single_time;
    idx_total = idx_total + 1; % Total iterations
    
    % Progress Display (Updates every 1000 iterations)
    if idx == 1000*idk
        Avg_Iteration_Time = 1000*Iteration_Time/idx;
        idk = idk + 1;
        disp(['Sources: ', num2str(source_idx)]);
        disp(['Iteration: ', num2str(idx), '  /  ', num2str(idx_total)]);
        disp(['1000 Iteration Time: ', num2str(Iteration_Time_1000), ' s']);
        disp(['Average 1000 Iteration Time: ', num2str(Avg_Iteration_Time), ' s']);
        if Iteration_Time < 60
            disp(['Running Model Time: ', num2str(Iteration_Time), ' s']);
        elseif Iteration_Time >= 60 && Iteration_Time < 3600
            disp(['Current Running Iteration Time: ', num2str(Iteration_Time/60), ' min']);
        else
            hour_single = floor(Iteration_Time/3600);
            disp(['Total Running Time: ', num2str(hour_single), ' hr, ',...
                num2str((Iteration_Time - (3600*hour_single))/60), ' min']);
        end
        if Total_Time < 3600
            disp(['Total Running Time: ', num2str(Total_Time/60), ' min']);
        else
            hour = floor(Total_Time/3600);
            disp(['Total Running Time: ', num2str(hour), ' hr, ',...
                num2str((Total_Time - (3600*hour))/60), ' min']);
        end
        disp(['# of Times Likelihood Increases: ',num2str(L_up_idx)]); %num2str(Loop_Sum)]);
        disp(['Similar Parameter Convergence:   ', num2str(Parameter_Break_VS_Likelihood_Break),...
            '     Amp: ', num2str(Sort_Amp_Std_Avg), '     Phi: ',...
            num2str(Sort_Phi_Std_Avg), '      Theta: ',...
            num2str(Sort_Theta_Std_Avg)  ]);
        disp(['# of Alternative Samples Used:   ', num2str(count_stuck)]);
        disp(['# of Random Sliced Samples Used: ', num2str(Random_Gen_Slice_Samp_idx)]);
        disp(['# of Direct Sliced Samples Used: ', num2str(Direct_Slice_Sample_idx)]);
        disp(['Stuck: ', num2str(stuck)]);
        % Likelihood curve plateau / Likelihood range of exploration table
        disp(['Threshold Check: ', num2str(Check_Threshold), ' / ', num2str(L_Dif_Initial_Table)]);
        % # of times likelihood does not meet convergence criteria / Exploration table convergence (1/0)
        disp(['Break: ', num2str(Loop_Sum), ' / ', num2str(Use_Other_Sampling_Methods)]);
        Iteration_Time_1000 = 0;
    end % End if display
    
% ----------------------------- Progress Plot -----------------------------
    % Plot likelhood as it develops along with current likelihood of
    % exploration table.
    if idx == Likelihood_Update*Run_plot_idx
        XAXIS = size(New_Table,1) + (1:size(Exploration_Table,1))-1;
        if Check_Threshold >= 5*Converge_Limit || L_Dif_Initial_Table >= 5*Exploration_Table_Threshold
            figure(102);
            plot(New_Table(:,1), 'r', 'linewidth', 2);
            hold on
            plot(XAXIS, sort(Exploration_Table(:,1)), 'b', 'linewidth', 2);
            hold off
            legend('New Table', 'Initial Table', 'location', 'northwest');
            title('Likelihood');
            xlabel('Iterations');
            ylabel('Likelihood');
        else
            figure(102);
            plot(New_Table(:,1), 'g', 'linewidth', 2);
            hold on
            plot(XAXIS, sort(Exploration_Table(:,1)), 'b', 'linewidth', 2);
            hold off
            legend('New Table', 'Initial Table', 'location', 'northwest');
            title('Likelihood');
            xlabel('Iterations');
            ylabel('Likelihood');
        end % End red vs green plot
        Run_plot_idx = Run_plot_idx + 1;
    end % End if plot
        
end % End While loop ------------------------------------------------------


%% CALCULATE EVIDENCE
% ------------------------- Evidence Calculations -------------------------
% Sort exploration table (lowest likelihood to highest)
Sorted_Table = sortrows(Exploration_Table,1);

% Append exploration table to the end of new table.
% Sorted from least to greatest by likelihood
% (Separated to keep all tables)
Final_Table = [New_Table; Sorted_Table]; 
switch source_idx
    case 1
        Final_Table_1 = Final_Table;
    case 2
        Final_Table_2 = Final_Table;
    case 3
        Final_Table_3 = Final_Table;
    case 4
        Final_Table_4 = Final_Table;
    case 5
        Final_Table_5 = Final_Table;
end
        

% For full likelihood curve (including the original table w/ replaced values)
Total_Likelihood = [New_Table(:,1); Sorted_Table(:,1)];

% Shift & unlog to prep for evidence calc
Final_L_log_noshift = New_Table(:,1); % Isolate likelihood

% Shift and unlog for evidence calcs
Shift_Final_L(source_idx) = max(Final_L_log_noshift); % Determine amount to shift lielihood
Sub_Shift_Final_L = Final_L_log_noshift - Shift_Final_L(source_idx);
Final_L = 10.^(Sub_Shift_Final_L); % Shift likelihood & unlog

% Initialize Bayesian iterative values
N = length(Final_L);
Z_sum = 0;

% Calculate Evidence (Eq. 20 in 'Bayesian Inference by Nested Sampling' - Tomislav Jasa and Ning Xiang)
for Z_idx = 1:N
    mu_E(Z_idx,source_idx) = exp(-Z_idx/N) - exp(-(Z_idx+1)/N);
    Like(Z_idx,source_idx) = Final_L(Z_idx);
    Z_nosum = Final_L(Z_idx) * mu_E(Z_idx,source_idx);
    Z_sum = Z_sum + Z_nosum;
end

% Finish log evidence calc
Z(source_idx) = Z_sum;
Evidence_log_noshift(source_idx) = log10(Z_sum);
Evidence_log(source_idx) = (Evidence_log_noshift(source_idx) + Shift_Final_L(source_idx));

% Calculate Bayes' Factor in log terms & decibans
if source_idx > 1
    B_i_j(source_idx) = Evidence_log(source_idx-1)/Evidence_log(source_idx);
    deciban(source_idx) = 10*log10(B_i_j(source_idx));
else
    B_i_j(source_idx) = Background_Noise_Evidence/Evidence_log(source_idx);
    deciban(source_idx) = 10*log10(B_i_j(source_idx));
end


%% DISPLAY SOURCE DoA ATTRIBUTES
% ------------------------ DoA Attributes Display -------------------------
% Create amplitude & direction matrix for all sources
par_idx_size = size(UI.amplitudes,2);
Amps(source_idx,:) = [UI.amplitudes zeros(1,num_sources_guess-par_idx_size)];
Phis(source_idx,:) = [UI.S_angles_az zeros(1,num_sources_guess-par_idx_size)];
Thetas(source_idx,:) = [UI.S_angles_el zeros(1,num_sources_guess-par_idx_size)];

disp(['Max Likelihood: ', num2str(max(Likelihood))]);
disp(['Max New Likelihood: ', num2str(max(New_Likelihood))]);
disp(' ');
disp(['Mod Amplitude: ', num2str(UI.amplitudes)]);
disp(['Mod Phi: ', num2str(UI.S_angles_az)]);
disp(['Mod Theta: ', num2str(UI.S_angles_el)]);

%% PLOT CURRENT EVIDENCE
% ---------------------------- Evidence Plots -----------------------------
% Plot Highest Likelihood Model
[~,Max_idx_L] = max(Final_Table(:,1));

% Reassign to rerun final calcs
UI.amplitudes = Final_Table(Max_idx_L,2:source_idx+1); % assign amplitude to rerun calcs
UI.S_angles_az = Final_Table(Max_idx_L,source_idx+2:2*source_idx+1); % assign source phi to rerun calcs
UI.S_angles_el = Final_Table(Max_idx_L,2*source_idx+2:3*source_idx+1); % assign source theta to rerun calcs

UI.plotflag = 1; % plot outputs
[Most_Likely_Model] = Beam_Model_Iterate_4_7_18(UI); % Beamform (create model)
title(strcat(num2str(source_idx), ' Source Model' ));
savefig(strcat(Save_Location,'\Beam', num2str(source_idx),'.fig'));

% Plot Final Likelihood
Like_XAXIS = size(New_Table,1) + (1:size(Sorted_Table,1))-1;
figure;
plot(New_Table(:,1), 'g', 'linewidth', 2);
hold on
plot(Like_XAXIS, Sorted_Table(:,1), 'r', 'linewidth', 2);
hold off
legend('Total Likelihood Curve', 'Exploration Table', 'location', 'northwest');
title(strcat(num2str(source_idx), ' Source Model Likelihood' ));
xlabel('Iterations');
ylabel('Likelihood');
savefig(strcat(Save_Location,'\L_i_', num2str(source_idx),'.fig'));

% Plot current log Evidence
figure(103);
bar(Evidence_log(1:source_idx))
title('Current log Evidence');
xlabel('Number of Sources');
ylabel('log Evidence');

% Plot current Evidence (in decibans)
figure(104);
bar(deciban(1:source_idx))
title('Current Deciban');
xlabel('Number of Sources');
ylabel('Decibans');


% =========================================================================
end % --------------------- End for Source Iterations ---------------------
% =========================================================================

%% FINAL PLOTS & DISPLAY
% ------------------------------ Final Plots ------------------------------
% Plot mu & Likliehood for each iteration
figure; hold on
for plot_idx = 1:num_sources_guess
        plot(1000*mu_E(:,plot_idx), 'linewidth', 1);
        plot(Like(:,plot_idx), 'linewidth', 1);
end
title('\mu(E(k)) & Likelihood');
legend('\mu_1', 'L_1', '\mu_2', 'L_2', '\mu_3', 'L_3', '\mu_4', 'L_4', ...
       '\mu_5', 'L_5', 'location', 'best');
savefig(strcat(Save_Location,'\mu(E(k)) & Likelihood.fig'));

% Plot Final Evidence
figure;
bar(Evidence_log)
title('log Evidence');
xlabel('Number of Sources');
ylabel('log Evidence');
savefig(strcat(Save_Location,'\log Evidence.fig'));

% Plot Final Decibans
figure;
bar(deciban)
title('Deciban Model Evidence');
xlabel('Number of Sources');
ylabel('Decibans');
savefig(strcat(Save_Location,'\Decibans.fig'));

% ----------------------------- Final Display -----------------------------
% Display parameters for each iteration
disp(' ');
disp(' ');
for a = 1:num_sources_guess
    disp([num2str(a),' -------------------------------------------------------------------------']);
    disp(['Mod Amplitude: ', num2str(Amps(a,:))]);
    disp(['Mod Phi: ', num2str(Phis(a,:))]);
    disp(['Mod Theta: ', num2str(Thetas(a,:))]);
    disp([num2str(a),' -------------------------------------------------------------------------']);
    disp(' ');
end

% Compile source attributes from lowest source number to highest for later viewing
Source_Attributes = [Amps Phis Thetas];

% Display total running time
if Total_Time < 60
    disp(['Total Running Time: ', num2str(Total_Time), ' s']);
elseif Total_Time >= 60 && Total_Time < 3600
    disp(['Total Running Time: ', num2str(Total_Time/60), ' min']);
else
    hour_single = floor(Total_Time/3600);
    disp(['Total Running Time: ', num2str(hour_single), ' hr, ',...
        num2str((Total_Time - (3600*hour_single))/60), ' min']);
end

% Save final workspace
save(strcat(Save_Location,'\',Load_Filename));
