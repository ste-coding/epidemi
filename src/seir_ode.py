# Modelo SEIR de campo médio para o sarampo, resolvido por RK4 escrito à mão.
# Limite de mistura homogênea -- Keeling & Rohani, Cap. 1-2.

from __future__ import annotations
import numpy as np

# parâmetros do sarampo (por dia): R0=15, período infeccioso 5 d, período latente 8 d
MEASLES = {"R0": 15.0, "gamma": 1.0 / 5.0, "sigma": 1.0 / 8.0}
MEASLES["beta"] = MEASLES["R0"] * MEASLES["gamma"]   # = 3,0 / dia


def seir_rhs(state, beta, sigma, gamma, N):
    S, E, I, R = state
    new_infections = beta * S * I / N
    dS = -new_infections
    dE = new_infections - sigma * E
    dI = sigma * E - gamma * I
    dR = gamma * I
    return np.array([dS, dE, dI, dR])


def rk4_step(state, dt, beta, sigma, gamma, N):
    k1 = seir_rhs(state, beta, sigma, gamma, N)
    k2 = seir_rhs(state + 0.5 * dt * k1, beta, sigma, gamma, N)
    k3 = seir_rhs(state + 0.5 * dt * k2, beta, sigma, gamma, N)
    k4 = seir_rhs(state + dt * k3, beta, sigma, gamma, N)
    return state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def solve_seir(N=10_000, I0=10, E0=0, days=120, dt=0.1, params=MEASLES):
    beta, sigma, gamma = params["beta"], params["sigma"], params["gamma"]
    state = np.array([N - I0 - E0, E0, I0, 0.0], dtype=float)
    n_steps = int(days / dt) + 1
    t = np.linspace(0, days, n_steps)
    Y = np.empty((n_steps, 4))
    Y[0] = state
    for k in range(1, n_steps):
        state = rk4_step(state, dt, beta, sigma, gamma, N)
        Y[k] = state
    return t, Y / N   # colunas: S, E, I, R (frações de N)


def final_size_and_peak(t, Y):
    I = Y[:, 2]
    attack_rate = Y[-1, 3]          # fração que já foi infectada
    peak_prev = I.max()
    peak_day = t[int(I.argmax())]
    return attack_rate, peak_prev, peak_day


if __name__ == "__main__":
    t, Y = solve_seir()
    ar, peak, day = final_size_and_peak(t, Y)
    print(f"R0 = {MEASLES['R0']:.1f}")
    print(f"beta,sigma,gamma = {MEASLES['beta']:.3f}, {MEASLES['sigma']:.3f}, {MEASLES['gamma']:.3f}")
    print(f"Attack rate = {ar*100:.1f}%")
    print(f"Peak prevalence = {peak*100:.2f}% on day {day:.1f}")
