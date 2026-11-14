import sympy as sp
v, mu, t = sp.symbols('v mu t', positive=True)


#D = sp.Function('D')(v, mu)

d = [sp.Function(f'd{k}')(v) for k in range(7)]
# D = d0(v) + d1(v) mu + ... + d6(v) mu^6
D = sum(d[k]*mu**k for k in range(7))

alpha = sp.Function('alpha')(v)
beta = sp.Function('beta')(v)
gamma = sp.Function('gamma')(v)
f0, f1, f2 = [sp.Function(n)(v, t) for n in ['f0', 'f1', 'f2']]
P = [sp.Integer(1), mu, (3*mu**2-1)/2]
# def L(h): return sp.diff(h, v)-mu/v*sp.diff(h, mu)
f = f0*P[0] + f1*P[1] + f2*P[2]


C = -(1/(v**2)) * sp.diff((v**2)*alpha*f, v) + \
    (1/(2*(v**2)))*sp.diff(sp.diff((v**2)*beta*f, v), v) + \
    (gamma/(4*(v**2)))*sp.diff((1-mu**2)*sp.diff(f, mu), mu)


G = (1-mu**2)*D*(sp.diff(f, v) - (mu/v)*sp.diff(f, mu))

Q = (1/(v**2)) * sp.diff((v**2) * G, v) - (1/v) * sp.diff(mu*G, mu)


df0_dt = (sp.Rational(1, 2)) * sp.integrate(P[0]*(C + Q), (mu, -1, 1))
df1_dt = (sp.Rational(3, 2)) * sp.integrate(P[1]*(C + Q), (mu, -1, 1))
df2_dt = (sp.Rational(5, 2)) * sp.integrate(P[2]*(C + Q), (mu, -1, 1))


df0_dt_c = sp.collect(sp.expand(df0_dt), mu)
df1_dt_c = sp.collect(sp.expand(df1_dt), mu)
df2_dt_c = sp.collect(sp.expand(df2_dt), mu)

print(f"df0_dt_c=")
sp.pprint(df0_dt_c)
print(sp.latex(df0_dt_c, mode="equation"))

print(f"df1_dt_c=")
sp.pprint(df1_dt_c)
print(sp.latex(df1_dt_c, mode="equation"))

print(f"df2_dt_c=")
sp.pprint(df2_dt_c)
print(sp.latex(df2_dt_c, mode="equation"))

"""
Cc = sp.collect(sp.expand(C), mu)
sp.pprint(Cc)
print(sp.latex(Cc, mode="equation"))

Qc = sp.collect(sp.expand(Q), mu)
sp.pprint(Qc)
print(sp.latex(Qc, mode="equation"))
"""
