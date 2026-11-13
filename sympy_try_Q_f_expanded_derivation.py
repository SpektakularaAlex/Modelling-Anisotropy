import sympy as sp
v, mu, t = sp.symbols('v mu t', positive=True)
D = sp.Function('D')(v, mu)
A, B, C = [sp.Function(n)(v, t) for n in ['A', 'B', 'C']]
P = [sp.Integer(1), mu, (3*mu**2-1)/2]
# def L(h): return sp.diff(h, v)-mu/v*sp.diff(h, mu)
f = A*P[0] + B*P[1] + C*P[2]

G = (1-mu**2)*D*(sp.diff(f, v) - (mu/v)*sp.diff(f, mu))

Q = (1/(v**2)) * sp.diff((v**2) * G, v) - (1/v) * sp.diff(mu*G, mu)

# sp.pprint(f"Q={Q}")

Qc = sp.collect(sp.expand(Q), mu)
sp.pprint(Qc)          # 2D pretty print in the terminal
print(sp.latex(Qc, mode="equation"))    # paste into Overleaf
