# Number of simulations 
n_sim = 40

# Random seed
seed = 23

# Firm innovation
F = 3
N0 = (2,2,2)
sigma2_eps = 0.0
theta = 1e-6
gamma = 0.0
ci = 0.0135

mu_ing = 0.0
mu_intact = -0.05
sigma2_ing = 0.01
sigma2_intact = 0.01

# Demand side
alpha = 0.5
cf = 0.1
m = 1.0
eta = 1.1

# Profit grid
q_gap = 2.0
nq = 21
timelimit_knitro = 30

muy_gap = 1.0
sigma2y_max = 0.2
nmuy = 7
nsigma2y = 6

# solve_approximation2
q_bound = 500.0
mu_bound1 = -0.25
mu_bound2 = -0.1
sigma_bound = 0.1
n_ing_cutoff = 200
# q_bound = 5.0
# mu_bound1 = -1.0
# mu_bound2 = -2.0
# sigma_bound = 1.0
# n_ing_cutoff = 2

# Integral
cutoff = 6
ngrid_int = 41

# Gurobi
timelimit = 600
mux_lb = -0.5
N_max = 8
threshold_sigma = 1e-5
bad_ing_bound = -0.05
integer_bound = 1e-5

# Simulation
T = 101

path = "simulation\parameters.txt"

with open(path, "w") as file:
    file.write("# Number of simulations \n")
    file.write(f"n_sim = {n_sim}\n\n")

    file.write("# Random seed\n")
    file.write(f"seed = {seed}\n\n")

    file.write("# Firm innovation\n")
    file.write(f"F = {F}\n")
    file.write(f"N0 = {N0}\n")
    file.write(f"sigma2_eps = {sigma2_eps}\n")
    file.write(f"theta = {theta}\n")
    file.write(f"gamma = {gamma}\n")
    file.write(f"ci = {ci}\n\n")
    file.write(f"mu_ing = {mu_ing}\n")
    file.write(f"mu_intact = {mu_intact}\n")
    file.write(f"sigma2_ing = {sigma2_ing}\n")
    file.write(f"sigma2_intact = {sigma2_intact}\n\n")

    file.write("# Demand side\n")
    file.write(f"alpha = {alpha}\n")
    file.write(f"cf = {cf}\n")
    file.write(f"m = {m}\n")
    file.write(f"eta = {eta}\n\n")

    file.write("# Profit grid\n")
    file.write(f"q_gap = {q_gap}\n")
    file.write(f"nq = {nq}\n")
    file.write(f"timelimit_knitro = {timelimit_knitro}\n\n")
    file.write(f"muy_gap = {muy_gap}\n")
    file.write(f"sigma2y_max = {sigma2y_max}\n")
    file.write(f"nmuy = {nmuy}\n")
    file.write(f"nsigma2y = {nsigma2y}\n")

    file.write("# solve_approximation2\n")
    file.write(f"q_bound = {q_bound}\n")
    file.write(f"mu_bound1 = {mu_bound1}\n")
    file.write(f"mu_bound2 = {mu_bound2}\n")
    file.write(f"sigma_bound = {sigma_bound}\n")
    file.write(f"n_ing_cutoff = {n_ing_cutoff}\n\n")

    file.write("# Integral\n")
    file.write(f"cutoff = {cutoff}\n")
    file.write(f"ngrid_int = {ngrid_int}\n\n")

    file.write("# Gurobi\n")
    file.write(f"timelimit = {timelimit}\n")
    file.write(f"mux_lb = {mux_lb}\n")
    file.write(f"N_max = {N_max}\n")
    file.write(f"threshold_sigma = {threshold_sigma}\n")
    file.write(f"bad_ing_bound = {bad_ing_bound}\n")
    file.write(f"integer_bound = {integer_bound}\n\n")

    file.write("# Simulation\n")
    file.write(f"T = {T}\n")
