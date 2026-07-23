# Dimensão de similaridade / fractal -- a ponte entre os dois cursos.
# Schiff, Cap. 1: D_s = log(N) / log(1/r). Curva de Koch: log4/log3 = 1,2619.
# Padrões epidêmicos reais são estatisticamente autossimilares, então estimamos a
# dimensão por contagem de caixas (igual a D_s para conjuntos exatamente autossimilares).

from __future__ import annotations
import numpy as np


def box_count(binary, box_sizes=None):
    # dimensão por contagem de caixas de um padrão booleano 2D
    binary = np.asarray(binary) > 0
    if binary.sum() == 0:
        return np.nan, np.array([]), np.array([])

    n = min(binary.shape)
    if box_sizes is None:
        max_p = int(np.floor(np.log2(n // 2))) if n >= 4 else 1
        box_sizes = [2 ** p for p in range(max_p + 1)]
    box_sizes = np.array(sorted(set(int(b) for b in box_sizes if b >= 1)))

    counts = []
    for b in box_sizes:
        pad_r = (-binary.shape[0]) % b
        pad_c = (-binary.shape[1]) % b
        padded = np.pad(binary, ((0, pad_r), (0, pad_c)), mode="constant")
        blocks = padded.reshape(padded.shape[0] // b, b, padded.shape[1] // b, b)
        counts.append(int(blocks.any(axis=(1, 3)).sum()))
    counts = np.array(counts)

    mask = counts > 0
    x = np.log(1.0 / box_sizes[mask])
    y = np.log(counts[mask])
    if len(x) < 2:
        return np.nan, box_sizes, counts
    dim = np.polyfit(x, y, 1)[0]   # inclinação = dimensão por contagem de caixas
    return float(dim), box_sizes, counts


def infection_dimension_over_time(frames, states=(1, 2)):
    # dimensão do conjunto ativo (E e I por padrão) em cada quadro
    out = {}
    for k, g in frames.items():
        out[k], _, _ = box_count(np.isin(g, states))
    return out


def koch_curve(order=6, length=1.0):
    p = np.array([[0.0, 0.0], [length, 0.0]])
    ang = -np.pi / 3
    rot = np.array([[np.cos(ang), -np.sin(ang)],
                    [np.sin(ang),  np.cos(ang)]])
    for _ in range(order):
        new = [p[0]]
        for i in range(len(p) - 1):
            a, b = p[i], p[i + 1]
            d = (b - a) / 3.0
            p1, p3 = a + d, a + 2 * d
            p2 = p1 + rot @ d           # ápice do triângulo equilátero
            new += [p1, p2, p3, b]
        p = np.array(new)
    return p[:, 0], p[:, 1]


def koch_similarity_dimension():
    return np.log(4) / np.log(3)


def rasterize_curve(x, y, grid=1024):
    # rasteriza a polilinha em uma grade booleana para a contagem de caixas
    img = np.zeros((grid, grid), dtype=bool)
    xn = (x - x.min()) / (x.max() - x.min() + 1e-12)
    yn = (y - y.min()) / (y.max() - y.min() + 1e-12)
    for i in range(len(xn) - 1):
        steps = max(2, int(grid * np.hypot(xn[i+1]-xn[i], yn[i+1]-yn[i])) + 1)
        tt = np.linspace(0, 1, steps)
        xs = ((xn[i] + tt * (xn[i+1] - xn[i])) * (grid - 1)).astype(int)
        ys = ((yn[i] + tt * (yn[i+1] - yn[i])) * (grid - 1)).astype(int)
        img[ys, xs] = True
    return img


if __name__ == "__main__":
    print(f"Koch exact D_s = log4/log3 = {koch_similarity_dimension():.5f}")
    x, y = koch_curve(order=6)
    d_box, _, _ = box_count(rasterize_curve(x, y, grid=1024))
    print(f"Koch box-counting estimate D_box = {d_box:.4f}")
