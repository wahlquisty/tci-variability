# Functions for TCI variability paper

struct PatientData # Patient covariates
    age::Float64
    bwt::Float64
    hgt::Float64
    opioid::Bool
    gdr::Int64
    av::Int64
end

struct PKPDparams # PKPD parameters
    V1::Float64
    V2::Float64
    V3::Float64
    CL::Float64
    Q2::Float64
    Q3::Float64
    Ce50::Float64
    gamma::Float64
    ke0::Float64
end


chi2cdf(n, d) = cdf(Chisq(n), d^2)  # Chi-squared cumulative distribution function

function eleveld_model(patient::PatientData; PKvariability=false, PDvariability=false, α=0.5)
    # Computes PK and PD parameters of the Eleveld covariate model for propofol.
    # Inputs: AGE (y), WGT (kg), HGT (cm), GDR (male 1, female 2).
    # add_opiates (false, no opiates. true opiates added)
    # A1V2 arterial or venous measurement (1 arterial, 2 venous)
    # PKvariability (if true, include PK variability)
    # PDvariability (if true, include PD variability)
    # α (quantile from which we draw our samples from the multivariate normal distribution)

    AGE = patient.age
    WGT = patient.bwt 
    HGT = patient.hgt
    GDR = patient.gdr
    add_opiates = patient.opioid
    A1V2 = patient.av

    # reference person (male, 35 yo, 170 cm, 70 kg)
    AGE_ref = 35
    WGT_ref = 70
    PMA_ref = (40 + AGE_ref * 52) / 52 # not born prematurely and now 35 yo
    BMI_ref = 24.2
    GDR_ref = 1 #1 male, 2 female

    BMI = WGT / ((HGT / 100)^2)
    PMA = (40 + AGE * 52) / 52 # Postmenstrual age

    # (reference) PK parameters
    theta_pk = [
        6.28        # V1ref [l]
        25.5        # V2ref [l]
        273         # V3ref [l]
        1.79        # Clref [l/min]
        1.83        # fel värde enlt corrigendum Q2ref [l/min]
        1.11        # Q3ref [l/min]
        0.191       # Typical residual error
        42.3        # CL maturation E50
        9.06        # CL maturation slope [weeks]
        -0.0156     # Smaller V2 with age
        -0.00286    # Lower CL with age
        33.6        # Weight for 50 # of maximal V1 [kg]
        -0.0138     # Smaller V3 with age
        68.3        # Maturation of Q3 [weeks]
        2.10        # CLref (female) [l/min]
        1.30        # Higher Q2 for maturation of Q3
        1.42        # V1 venous samples (children)
        0.68        # Higer Q2 venous samples
    ]

    # (reference) PD model
    theta_pd = [
        3.08       # Ce50 [ug/ml]
        0.146      # ke0 for arterial samples [1/min]
        93.0       # Baseline BIS value
        1.47       # PD sigmoid slope (Ce > Ce50)
        8.03       # Residual error [BIS]
        0.0517     # Increase in delay with age
        -0.00635   # Decrease in Ce50 with age
        1.24       # ke0 for venous samples [1/min]
        1.89]      # PD sigmoid slope (Ce < Ce50)


    # Inter-individual variability
    if PKvariability || PDvariability # If there is any variability in the PK or PD parameters


        n_eta_pk = 0
        n_eta_pd = 0
        sigma_pk = []
        sigma_pd = []

        if PKvariability # If there is any variability in thePK parameters
            # variability for PK
            variance_pk = [
                0.610
                0.565
                0.597
                0.265
                0.346
                0.209]
            # 0.463] # Residual error, not used

            sigma_pk = sqrt.(variance_pk) # Standard deviations
            n_eta_pk = length(variance_pk)
        end

        if PDvariability
            variance_pd = [
                0.242      # Ce50
                0.702]      # ke0
            # 0.230]     # Residual error

            sigma_pd = sqrt.(variance_pd)
            n_eta_pd = length(variance_pd)
        end

        n_eta = n_eta_pk + n_eta_pd
        etas = zeros(n_eta) # final etas
        sigmas = [sigma_pk; sigma_pd]

        # Draw eta from a normal distribution with mean 0 and standard deviation η_i. If not within α quantile, discard and draw new parameter values

        Q = zeros(n_eta, n_eta)
        for j = 1:n_eta # Q = [σ1^2 0; 0 σ2^2]
            Q[j, j] = sigmas[j]^2
        end
        # Normal distribution where α is the probability for being within the set of patients closest to the mean
        distr = Vector{typeof(Normal(0, 1))}(undef, n_eta) # create distributions
        for j = 1:n_eta
            distr[j] = Normal(0, sigmas[j]) # Create normal distribution
        end
        # While loop
        x = zeros(n_eta)
        pkeep = 1.1 # the probability of this sample to keep, initialize at some value outside of likely probabilities
        while pkeep > α
            for j = 1:n_eta
                x[j] = rand(distr[j]) # Draw random values from the distributions
            end
            d = mahalanobis(x, zeros(n_eta), inv(Q))
            pkeep = chi2cdf(n_eta, d)
            # println(α)
        end
        etas = x
    end

    if PKvariability && PDvariability # Get final etas
        eta_pk = etas[1:n_eta_pk]
        eta_pd = etas[n_eta_pk+1:end]
    elseif PKvariability
        eta_pk = etas
        eta_pd = zeros(2)
    elseif PDvariability
        eta_pk = zeros(6)
        eta_pd = etas
    else
        eta_pk = zeros(6)
        eta_pd = zeros(2)
    end

    # helper functions
    f_ageing(x) = exp(x * (AGE - AGE_ref))
    f_sigmoid(x, E50, lambda) = (x^lambda) / (x^lambda + E50^lambda)

    function f_AlSallami(GDR, AGE, WGT, BMI)
        if GDR == 1 #male
            al = (0.88 + (1 - 0.88) / (1 + (AGE / 13.4)^(-12.7))) * (9270 * WGT) / (6680 + 216 * BMI)
        else # female
            al = (1.11 + (1 - 1.11) / (1 + (AGE / 7.1)^(-1.1))) * (9270 * WGT) / (8780 + 244 * BMI)
        end
        return al
    end

    f_central(x) = f_sigmoid(x, theta_pk[12], 1)
    f_CLMaturation = f_sigmoid(PMA * 52, theta_pk[8], theta_pk[9])
    f_CLMaturation_ref = f_sigmoid(PMA_ref * 52, theta_pk[8], theta_pk[9])
    f_Q3Maturation = f_sigmoid(AGE * 52 + 40, theta_pk[14], 1)
    f_Q3Maturation_ref = f_sigmoid(AGE_ref * 52 + 40, theta_pk[14], 1)
    function f_opiates(x)
        if add_opiates == true
            op = exp(x * AGE)
        else
            op = 1              # due to absence of opiates
        end
        return op
    end
    f_Al = f_AlSallami(GDR, AGE, WGT, BMI)
    f_Al_ref = f_AlSallami(GDR_ref, AGE_ref, WGT_ref, BMI_ref)

    # compartment volumes [L]
    if A1V2 == 1 # arterial
        V1 = theta_pk[1] * f_central(WGT) / f_central(WGT_ref) * exp(eta_pk[1])
    else
        V1 = theta_pk[1] * f_central(WGT) / f_central(WGT_ref) * exp(eta_pk[1]) * (1 + theta_pk[17] * (1 - f_central(WGT)))
    end

    V2 = theta_pk[2] * (WGT / WGT_ref) * f_ageing(theta_pk[10]) * exp(eta_pk[2])
    V2_ref = theta_pk[2]
    V3 = theta_pk[3] * (f_Al / f_Al_ref) * f_opiates(theta_pk[13]) * exp(eta_pk[3])
    V3_ref = theta_pk[3]

    # clearances [L/min]
    CL = ((2 - GDR) * theta_pk[4] + (GDR - 1) * theta_pk[15]) * ((WGT / WGT_ref)^(0.75)) * (f_CLMaturation / f_CLMaturation_ref) * f_opiates(theta_pk[11]) * exp(eta_pk[4])  # elimination clearance
    if A1V2 == 1 # arterial
        Q2 = theta_pk[5] * (V2 / V2_ref)^(0.75) * (1 + theta_pk[16] * (1 - f_Q3Maturation)) * exp(eta_pk[5])    # compartment clearance
    else
        Q2 = theta_pk[18] * theta_pk[5] * (V2 / V2_ref)^(0.75) * (1 + theta_pk[16] * (1 - f_Q3Maturation)) * exp(eta_pk[5]) #compartment clearance
    end
    Q3 = theta_pk[6] * (V3 / V3_ref)^(0.75) * (f_Q3Maturation / f_Q3Maturation_ref) * exp(eta_pk[6])      # compartment clearance


    if A1V2 == 1 # arterial measurement
        ke0 = theta_pd[2] * (WGT / WGT_ref)^(-0.25) * exp(eta_pd[1]) # [1/min]
    else # venous measurement
        ke0 = theta_pd[8] * (WGT / WGT_ref)^(-0.25) * exp(eta_pd[2]) # [1/min]
    end

    Ce50 = theta_pd[1] * f_ageing(theta_pd[7]) * exp(eta_pd[1])
    #Ce50 = Ce50*V1; # scale from conc to mass
    gamma = theta_pd[4] # Assume Ce <= Ce50
    # BIS_var = theta_pd[5] * exp(eta_pd[3])

    pkpdparams = PKPDparams(V1, V2, V3, CL, Q2, Q3, Ce50, gamma, ke0)
    # return V1, V2, V3, CL, Q2, Q3, Ce50, gamma, ke0
    return pkpdparams
end

function createlinearPKPDmodel(params::PKPDparams)
    V1 = params.V1
    V2 = params.V2
    V3 = params.V3
    CL = params.CL
    Q2 = params.Q2
    Q3 = params.Q3
    ke0 = params.ke0

    k10 = CL / V1  #[min^-1]
    k12 = Q2 / V1  #[min^-1]
    k13 = Q3 / V1  #[min^-1]
    k21 = Q2 / V2  #[min^-1]
    k31 = Q3 / V3  #[min^-1]

    # A matrix for pk model
    A_pk = zeros(3, 3)
    A_pk[1, 1] = -(k10 + k12 + k13)
    A_pk[1, 2] = k12
    A_pk[1, 3] = k13
    A_pk[2, 1] = k21
    A_pk[2, 2] = -(k21)
    A_pk[3, 1] = k31
    A_pk[3, 3] = -(k31)

    # B matrix for pk model
    B_pk = zeros(3)
    B_pk[1, 1] = 1 / V1

    # Define the effect-site PD and combine

    Te = 1 / ke0 # Effect site time constant [min] 90-100 sec

    # A and B matrices for pkpd-model
    A_pkpd = [A_pk zeros(3, 1); 1/Te zeros(1, 2) -1/Te]
    B_pkpd = [B_pk; 0]

    return A_pkpd, B_pkpd
end

function discretization(A_pkpd, B_pkpd,h)
    n = 4                    # Number of states
    M = expm([A_pkpd B_pkpd; zeros(1, n + 1)] * h)
    Phi = M[1:n, 1:n]
    Gamma = M[1:n, n+1]
    return n, Phi, Gamma
end

function getlimitrefconc(Ce50, gamma)
    BIS_baseline = 100 #93

    Ce_BIS40 = ((BIS_baseline - 40) / 40)^(1 / gamma) * Ce50
    Ce_BIS50 = ((BIS_baseline - 50) / 50)^(1 / gamma) * Ce50

    CeMax = Ce_BIS40 * ones(N)
    r = Ce_BIS50 * ones(N)

    return BIS_baseline, CeMax, r
end

eye(n) = Matrix(1.0I, n, n)

function QP(Phi, Gamma, CeMax, r, N, x0) # Set up the QP
    n = length(x0)
    Phij = eye(n) # Phi^j
    Phij1Gamma = zeros(N) # Row j will hold Phi^j_1*Gamma
    Phij4Gamma = zeros(N) # Row j will hold Phi^j_4*Gamma
    E1 = zeros(N, n)
    E4 = zeros(N, n)
    F1 = zeros(N, N)
    F4 = zeros(N, N)
    for j = 1:N
        Phij1Gamma[j] = Phij[1, :]' * Gamma
        Phij4Gamma[j] = Phij[4, :]' * Gamma
        for i = 1:j
            F1[j, i] = Phij1Gamma[j-i+1]
            F4[j, i] = Phij4Gamma[j-i+1]
        end
        Phij = Phij * Phi
        E1[j, :] = Phij[1, :]
        E4[j, :] = Phij[4, :]
    end

    D = eye(N)
    # D[end, end] = alpha
    # Speed-up is possible here:
    # 1. H is symmetric, so only need to compute half
    # 2. D and F1 are both sparse
    DF4 = D * F4                     # Compute only once
    #H=F1'*DF1                   # Compute only once
    H = sparse(0.5 * (F4' * DF4 + DF4' * F4))      # Same but numerically more robust

    f0T = E4' * DF4                  # Compute only once
    f1T = r' * DF4                   # Recompute in every iteration
    fT = x0' * f0T - f1T               # Recompute in every iteration
    f = fT'
    A = sparse([-eye(N); F4])               # Compute only once
    b = vec([zeros(N, 1); CeMax - E4 * x0])  # Recompute lower part in each iteration

    return H, f, A, b
end

function simulate(Phi, Gamma, u, h, N, x0)
    t = h .* (0:1:N-1) # Time vector
    n = length(x0)

    # Simulate the system
    xhats = zeros(N, n)
    xhats[1, :] = x0
    xhatprev = xhats[1, :]
    for i = 2:N
        xhatnew = Phi * xhatprev + Gamma * u[i-1]
        xhats[i, :] = xhatnew
        xhatprev = xhatnew
    end
    return t, xhats
end

computeBIS(x, BIS_baseline, Ce50, gamma) = BIS_baseline .* (Ce50 .^ gamma ./ (Ce50 .^ gamma .+ x .^ gamma))

function getrefpatient(h)
    age_ref = 35 # Patient age [years]
    bwt_ref = 70 # Patient weight [kg]
    hgt_ref = 170 # Patient height [cm]
    opioid_ref = false # Model not based on added remifentanil
    gdr_ref = 1 # Model based on male
    av_ref = 1  # Model based on venous sampling

    referencecovariates = PatientData(age_ref, bwt_ref, hgt_ref, opioid_ref, gdr_ref, av_ref)

    # Create reference patient, with all eta put to zero
    pkpdparams = eleveld_model(referencecovariates, PKvariability=false, PDvariability=false) # Get all parameters for volumes and clearances from the Eleveld model
    A_pkpd, B_pkpd = createlinearPKPDmodel(pkpdparams) # Get A and B matrices for the pkpd model
    n, Phi, Gamma = discretization(A_pkpd, B_pkpd, h) # ZOH discretization

    return referencecovariates, pkpdparams, n, Phi, Gamma
end

function TCIreference(pkpdparams, n, Phi, Gamma, N, h)
    
    BIS_baseline, CeMax, r = getlimitrefconc(pkpdparams.Ce50, pkpdparams.gamma)
    x0 = zeros(n) # Initial state
    H, f, A, b = QP(Phi, Gamma, CeMax, r, N, x0) # Set up the QP

    # Solve the QP
    prob = OSQP.Model()
    OSQP.setup!(prob; P=H, q=f, A=A, l=nothing, u=b,
        adaptive_rho=true # verkar inte göra någon skillnad
    )
    results = OSQP.solve!(prob)

    # Extract the solution
    u = results.x

    # Plot the BIS for the reference patient
    t, X = simulate(Phi, Gamma, u, h, N, x0)
    x4 = X[:, 4] # Fourth state, effect site concentration
    BIS = computeBIS(x4, BIS_baseline, pkpdparams.Ce50, pkpdparams.gamma) # Convert to BIS

    return t, u, BIS, N, r, CeMax, BIS_baseline
end

function plotBISu(t, u, BIS)
    p1 = plot(t, BIS, label="BIS", xlabel="Time (min)", ylabel="BIS") # Plot BIS for reference patient
    p2 = plot(t, u, label="Dose", xlabel="Time (min)", ylabel="Dose (mg/min)")
    plot(p1, p2, layout=(2, 1), legend=:topright)
end

# Use same values as reference patient, but include the inter-individual
#variability eta and simulate response for n_pat patients using the same u.
function samplesimulatepatients(h, covariates, BIS_baseline, n_pat_sample, α, PKvariability, PDvariability; seed=123456)
    Random.seed!(seed)
    X4_variability = zeros(length(t), n_pat_sample) # Save X4 values for all patients
    BIS_all = similar(X4_variability) # Save BIS values for all patients
    for j = 1:n_pat_sample
        pkpdparams_j = eleveld_model(covariates, PKvariability=PKvariability, PDvariability=PDvariability; α=α)
        A_pkpd, B_pkpd = createlinearPKPDmodel(pkpdparams_j)
        n, Phi, Gamma = discretization(A_pkpd, B_pkpd, h)
        x0 = zeros(n) # Initial state
        t, X4_variability = simulate(Phi, Gamma, u, h, N, x0)
        BIS_all[:, j] = computeBIS(X4_variability[:, 4], BIS_baseline, pkpdparams_j.Ce50, pkpdparams_j.gamma)
    end
    return BIS_all
end

function plotsampledpatients(n_pat, BIS_all, t, BIS_ref)     # Plot simulation results
    p = plot(xlabel="Time (min)", ylabel="BIS", ylims=(0, 100))
    for i = 1:n_pat
        plot!(p, t, BIS_all[:, i], label="", color=1)
    end
    plot!(p, t, BIS_ref, label="Ref", color=2, width=5)
    display(p)
end


function saveresultstofile(α, t, u, BIS_all, BIS_ref)
    fname = string(Int(α*100))

    # Save results to csv file
    t_final = Vector(t)
    y_final = BIS_all'
    data = [t_final'; y_final]'
    df = DataFrame(data, :auto)
    CSV.write("csv/tci-variability-" * fname * "percent.csv", df)

    # Reference patient data
    y_ref = BIS_ref
    data_ref = [t'; y_ref']'
    df_ref = DataFrame(data_ref, :auto)
    CSV.write("csv/bisref.csv", df_ref)

    # Save dose to csv file
    u_final = abs.(Vector(u))
    data_dose = [t_final'; u_final']'
    df_dose = DataFrame(data_dose, :auto)
    CSV.write("csv/dose.csv", df_dose)
end


function savehistdatatofile(BIS_all, α, n_pat_sample) # Save histogram data of steady state BIS values
    fname = string(Int(α*100))

    # Create hist data of final BIS values
    BISend = BIS_all[end, :]
    BISendint = Int64.(round.(BISend, digits=0))
    A = BISendint
    BISendcount = hcat([[i; count(==(i), A)] for i in unique(A)]...)' # count how many of each value

    B = BISendcount[sortperm(BISendcount[:, 1]), :]

    # Fill out with zeros for missing values
    for i = 1:100
        if i ∉ B[:, 1]
            # @show i
            B = vcat(B, [i 0])
        end
    end

    # Sort by BIS value
    B = B[sortperm(B[:, 1]), :]

    # Collect five BIS values together to get a smoother histogram
    C = reshape(B[:, 1], (5, 20))'
    # Choose the average BIS value for each group of five
    D = vec(sum(C, dims=2) ./ 5)
    C2 = reshape(B[:, 2], (5, 20))'
    # Compute sum of counts for each group of five
    D2 = vec(sum(C2, dims=2))

    plot(D, D2, label="", xlabel="BIS", ylabel="Count", legend=:none)

    # Save hist data of final BIS values
    df_hist = DataFrame(:BIS => D, :Count => D2)
    CSV.write("csv/bisend" * fname * "hist" * string(n_pat_sample) * ".csv", df_hist)
end

function countBISoutside(BIS_all) # Count number of patients with BIS > 60 or BIS < 40 at the end of the simulation
    BISend = BIS_all[end, :]
    BISendint = Int64.(round.(BISend, digits=0))
    A = BISendint
    BISendcount = hcat([[i; count(==(i), A)] for i in unique(A)]...)' # count how many of each value

    B = BISendcount[sortperm(BISendcount[:, 1]), :]

    # Count number of patients with BIS > 60 or BIS < 40 at the end of the simulation
    num60 = B[findall(B[:, 1] .> 60), 2] |> sum # Number of patients with BIS > 60
    num40 = B[findall(B[:, 1] .< 40), 2] |> sum # Number of patients with BIS < 40
    totaloutside = num60 + num40 # Total number of patients with BIS outside the interval [40, 60]

    return num60, num40, totaloutside
end