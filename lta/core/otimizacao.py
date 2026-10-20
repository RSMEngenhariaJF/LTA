"""Otimização sem derivadas (Nelder-Mead), para ajustes pequenos sem depender do SciPy."""
import numpy as np


def nelder_mead(func, x0, passo=0.2, tol=1e-12, max_iter=4000):
    """Minimiza func(x) a partir de x0 (vetor). Retorna (x, f(x))."""
    x0 = np.asarray(x0, dtype=float)
    n = len(x0)
    pts = [x0] + [x0 + passo * np.eye(n)[i] for i in range(n)]
    val = [func(p) for p in pts]
    for _ in range(max_iter):
        ordem = np.argsort(val)
        pts = [pts[i] for i in ordem]
        val = [val[i] for i in ordem]
        if abs(val[-1] - val[0]) < tol and np.max(np.abs(np.array(pts) - pts[0])) < 1e-10:
            break
        centro = np.mean(pts[:-1], axis=0)
        xr = centro + (centro - pts[-1])
        fr = func(xr)
        if fr < val[0]:
            xe = centro + 2 * (centro - pts[-1])
            fe = func(xe)
            pts[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            pts[-1], val[-1] = xr, fr
        else:
            xc = centro + 0.5 * (pts[-1] - centro)
            fc = func(xc)
            if fc < val[-1]:
                pts[-1], val[-1] = xc, fc
            else:
                pts = [pts[0]] + [pts[0] + 0.5 * (p - pts[0]) for p in pts[1:]]
                val = [val[0]] + [func(p) for p in pts[1:]]
    i = int(np.argmin(val))
    return pts[i], val[i]
