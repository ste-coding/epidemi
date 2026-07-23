# Autômato celular SEIR 2D probabilístico para o sarampo.
# Contraparte espacial de contato local da EDO SEIR de campo médio.
# Fuentes & Kuperman, Physica A 267 (1999); Schiff, Cap. 1.

from __future__ import annotations
import numpy as np
from scipy.ndimage import convolve

# códigos de estado
S, E, I, R = 0, 1, 2, 3

MOORE = np.array([[1, 1, 1],
                  [1, 0, 1],
                  [1, 1, 1]], dtype=int)

VON_NEUMANN = np.array([[0, 1, 0],
                        [1, 0, 1],
                        [0, 1, 0]], dtype=int)


class SEIR_CA:
    def __init__(self, L=200, tau=0.30, sigma=1/8, gamma=1/5,
                 neighbourhood="moore", boundary="constant", seed=0):
        self.L = L
        self.tau = tau
        self.sigma = sigma
        self.gamma = gamma
        self.kernel = MOORE if neighbourhood == "moore" else VON_NEUMANN
        self.z = int(self.kernel.sum())          # número de coordenação
        self.boundary = boundary                 # 'constant' (fechada) ou 'wrap'
        self.rng = np.random.default_rng(seed)

        # todos suscetíveis, exceto uma única semente infecciosa no centro
        self.grid = np.zeros((L, L), dtype=np.int8)
        self.grid[L // 2, L // 2] = I

    def local_R0_proxy(self):
        # proxy de reprodução em mistura homogênea; o R0 real no reticulado é menor
        return self.z * self.tau / self.gamma

    def counts(self):
        g = self.grid
        return ((g == S).sum(), (g == E).sum(), (g == I).sum(), (g == R).sum())

    def step(self):
        g = self.grid
        infected = (g == I).astype(np.int8)
        n_inf = convolve(infected, self.kernel, mode=self.boundary, cval=0)

        new = g.copy()
        # S -> E conforme o número de vizinhos infecciosos
        r = self.rng.random(g.shape)
        p_inf = 1.0 - (1.0 - self.tau) ** n_inf
        new[(g == S) & (r < p_inf)] = E
        # E -> I
        r2 = self.rng.random(g.shape)
        new[(g == E) & (r2 < self.sigma)] = I
        # I -> R
        r3 = self.rng.random(g.shape)
        new[(g == I) & (r3 < self.gamma)] = R
        self.grid = new

    def run(self, steps=200, record_frames=None):
        record_frames = record_frames or []
        Ncells = self.L * self.L
        ts = np.empty((steps + 1, 4))
        ts[0] = np.array(self.counts()) / Ncells
        frames = {}
        if 0 in record_frames:
            frames[0] = self.grid.copy()
        for k in range(1, steps + 1):
            self.step()
            ts[k] = np.array(self.counts()) / Ncells
            if k in record_frames:
                frames[k] = self.grid.copy()
            # para quando a epidemia se extingue
            if ts[k, 1] == 0 and ts[k, 2] == 0:
                ts = ts[:k + 1]
                break
        return ts, frames


if __name__ == "__main__":
    ca = SEIR_CA(L=200, tau=0.30, seed=1)
    print(f"coordination z = {ca.z}")
    print(f"local R0 proxy (z*tau/gamma) = {ca.local_R0_proxy():.1f}")
    ts, _ = ca.run(steps=200)
    print(f"CA peak prevalence = {ts[:,2].max()*100:.2f}% on step {int(ts[:,2].argmax())}")
    print(f"CA attack rate = {ts[-1,3]*100:.1f}%")
