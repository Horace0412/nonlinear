import numpy as np

import time
import random

import config

from gurobipy import Model, GRB
import gurobipy as gp

from knitro import *

from scipy.io import savemat
from scipy.io import loadmat

def mu_sigma_ridbad(NN, iing, pp, mu, sigma):
    if NN <= 6:
        return NN, iing, pp, mu, sigma, np.where(mu[:NN] >= -1e6)[0]

    good_ing = np.where(mu[:NN] >= config.bad_ing_bound)[0]
    bad_ing = set(np.where(mu[:NN] < config.bad_ing_bound)[0])
    
    NN_ridbad = NN - len(bad_ing)
    iing_ridbad = iing[good_ing]
    pp_ridbad = NN_ridbad + NN_ridbad*(NN_ridbad-1)//2

    bad_ing_m_index = []

    for i in range(NN):
        if i in bad_ing:
            bad_ing_m_index.append(i)

    offset = 0
    for i in range(NN):
        for j in range(i+1, NN):
            m_index = NN + offset
            if i in bad_ing or j in bad_ing:
                bad_ing_m_index.append(m_index)
            offset += 1

    mu_ridbad = np.delete(mu, bad_ing_m_index, axis=0)
    sigma_ridbad_temp = np.delete(sigma, bad_ing_m_index, axis=0)
    sigma_ridbad = np.delete(sigma_ridbad_temp, bad_ing_m_index, axis=1)

    return NN_ridbad, iing_ridbad, pp_ridbad, mu_ridbad, sigma_ridbad, good_ing










def build_variance_expr(N, p):
    # create a new Gurobi model
    model = gp.Model("extended_variable_formulation")

    # create p binary variables m
    m = model.addVars(p, vtype=GRB.BINARY, name="m")
    
    # enforce m[N + offset] = x_i*x_j for each pair (i<j)
    offset_map = {}
    offset = 0
    for i in range(N):
        for j in range(i+1, N):
            m_index = N + offset
            offset_map[(i,j)] = m_index
            # constraints for binary product
            model.addConstr(m[m_index] <= m[i], name=f"int_le_i[{i},{j}]")
            model.addConstr(m[m_index] <= m[j], name=f"int_le_j[{i},{j}]")
            model.addConstr(m[m_index] >= m[i] + m[j] - 1, name=f"int_ge[{i},{j}]")
            offset += 1

    return model, m, offset_map

def mu_exp_gen(p, m, mu):
    mu_exp = gp.LinExpr()
    for k in range(p):
        mu_exp.addTerms(mu[k], m[k])
    return mu_exp

def sigma2_exp_gen(model, p, m, sigma): 
    # this function computes mT*sigma*m, using Cholesky decomposition to speed up
    L = np.linalg.cholesky(sigma + 1e-6 * np.eye(sigma.shape[0]))
    z = model.addVars(p, lb = -100, vtype=GRB.CONTINUOUS, name="z")
    for r in range(p):
        model.addConstr(
            z[r] == gp.quicksum(L[j, r] * m[j] for j in range(p) if abs(L[j, r]) > config.threshold_sigma),
            name=f"z_eq_proj[{r}]"
        )
    sigma2_exp = gp.quicksum(z[r] * z[r] for r in range(p))

    return sigma2_exp

def gurobi_muSigma_sample(N, p, mu, sigma, q, timelimit, muy_grid, sigma2y_grid, profits_grid2):
    # create Gurobi model
    model, m, offset_map = build_variance_expr(N, p)
    
    mux_exp = mu_exp_gen(p, m, mu)
    sigma2x_exp = config.sigma2_eps + sigma2_exp_gen(model, p, m, sigma)

    # define mux
    mux = model.addVar(vtype=GRB.CONTINUOUS, lb = q + config.mux_lb, name="mux")
    model.addConstr(mux == mux_exp, name="con_mux")

    # define sigmax
    sigma2x = model.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name="sigma2x")
    model.addConstr(sigma2x == sigma2x_exp, name="con_sigmax")

    # define penalty term: theta * (sum x)^gamma
    ning_grid = list(range(N + 1))
    cost_grid = [config.theta * k**config.gamma for k in ning_grid]
    s = model.addVar(vtype=GRB.INTEGER, lb=0, ub=N, name="s")
    t = model.addVar(vtype=GRB.CONTINUOUS, name="t", lb=0.0)
    model.addConstr(s == gp.quicksum(m[i] for i in range(N)), name="sum_m_main_effects")
    model.addConstr(s >= 1, name="non_empty_set")
    model.addConstr(s <= config.N_max, name="optimality")
    model.addGenConstrPWL(s, t, ning_grid, cost_grid, name="pwl_power")

    # SOS2 constraints
    nx, ny = len(muy_grid), len(sigma2y_grid)
    lam_x = model.addVars(nx, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="lam_x")
    lam_y = model.addVars(ny, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="lam_y")

    model.addConstr(gp.quicksum(lam_x[i] for i in range(nx)) == 1)
    model.addConstr(gp.quicksum(lam_y[j] for j in range(ny)) == 1)

    model.addConstr(mux == gp.quicksum(lam_x[i]*muy_grid[i] for i in range(nx))) # enforce interpolation constraints (x,y interpolation)
    model.addConstr(sigma2x == gp.quicksum(lam_y[j]*sigma2y_grid[j] for j in range(ny)))
    
    model.addSOS(GRB.SOS_TYPE2, [lam_x[i] for i in range(nx)])
    model.addSOS(GRB.SOS_TYPE2, [lam_y[j] for j in range(ny)])

    # set the objective function (maximize)
    model.setObjective(gp.quicksum(lam_x[i] * lam_y[j] *profits_grid2[i,j] for i in range(nx) for j in range(ny)) - t, GRB.MAXIMIZE)
    
    # model.setParam("NonConvex", 2)
    model.setParam("OutputFlag", 0)
    # model.setParam("TimeLimit", timelimit)
    # model.setParam("Method", 0)
    # model.setParam("Cuts", 0)
    # model.setParam("PrePasses", 5)

    # solve the model
    model.optimize()

    # extract solution
    x_opt = np.array([m[i].X for i in range(N)])

    # print([(v.VarName, v.X) for v in model.getVars()])

    # model analysis
    # model.tune()
    # model.getTuneResult(0)

    return model.objVal, x_opt, model.status, model.MIPGap








for i_sim in range(2,3):
    print(f"Simulation No. {i_sim}")

    # set random seed to guarantee reproductibility
    np.random.seed(config.seed + i_sim)            # NumPy seed
    random.seed(config.seed + i_sim)               # random seed
    rng = np.random.default_rng(config.seed + i_sim)

    data = loadmat(fr"simulation\result_{i_sim}.mat")

    # read data
    turn = data["turn"].squeeze()
    N = data["N"].squeeze()

    N_f = data["N_f"].squeeze()
    ing_f = data["ing_f"].squeeze()
    p_f = data["p_f"].squeeze()
    mu_f = data["mu_f"].squeeze()
    sigma_f = data["sigma_f"].squeeze()
    q_beg_f = data["q_beg_f"].squeeze()
    q_end_f = data["q_end_f"].squeeze()

    newing_f = data["newing_f"].squeeze()
    newdrug_f = data["newdrug_f"].squeeze()
    x_f = data["x_f"].squeeze()
    x_n_f = data["x_n_f"].squeeze()
    v_f = data["v_f"].squeeze()
    x_ing_f = data["x_ing_f"].squeeze()
    x_ing_n_f = data["x_ing_n_f"].squeeze()
    v_ing_f = data["v_ing_f"].squeeze()
    y_f = data["y_f"].squeeze()
    eps = data["eps"].squeeze()

    price_f = data["price_f"].squeeze()
    share_f = data["share_f"].squeeze()
    profits_f = data["profits_f"].squeeze()

    knitro_status = data["knitro_status"].squeeze()
    knitro_status_comp = data["knitro_status_comp"].squeeze()
    sol_status = data["sol_status"].squeeze()
    sol_status_ing = data["sol_status_ing"].squeeze()
    MIPGap = data["MIPGap"].squeeze()
    MIPGap_ing = data["MIPGap_ing"].squeeze()
    time_taken = data["time_taken"].squeeze()
    total_time_taken = data["total_time_taken"].squeeze()

    alpha_truth = data["alpha_truth"].squeeze()




    t = config.T - 2
    N[t] = len(np.unique(np.concatenate([ing_f[t][0][0],ing_f[t][1][0],ing_f[t][2][0]])))
    active_firm = turn[t]

    # A. read grid
    data = loadmat(fr"simulation\output_data.mat")
    muy_grid = data["muy_grid"].squeeze()
    sigma2y_grid = data["sigma2y_grid"].squeeze()
    profits_grid2 = data["profits_grid2"].squeeze()


    # B. optimization
    NN, iing, pp, mu, sigma, q = N_f[t,active_firm].item(), ing_f[t][active_firm][0], p_f[t,active_firm].item(), mu_f[t][active_firm][0], sigma_f[t][active_firm], q_beg_f[t,active_firm]
    NN_ridbad, iing_ridbad, pp_ridbad, mu_ridbad, sigma_ridbad, good_ing = mu_sigma_ridbad(NN, iing, pp, mu, sigma)
    a = time.time()
    v, x_opt_temp, sol_status, MIPGap = gurobi_muSigma_sample(NN, pp, mu, sigma, q, config.timelimit, muy_grid, sigma2y_grid, profits_grid2)
    b = time.time()
    print(b-a)


    