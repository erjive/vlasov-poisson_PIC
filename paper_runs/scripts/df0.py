"""F0(Q,J) en Python, replicando src/distribution.f90 linea por linea."""
import numpy as np
L0 = 2.0
c  = 0.5*(L0 + np.sqrt(L0**2 + 4.0))

BIM = dict(b1=0.8, b2=0.3, Ja=0.10, sa=0.025, Jb=0.24, sb=0.035, wb=0.7)
SPI = dict(J0=0.15, s0=0.04, beta=50.0, sq=0.5)
KIN = dict(Jt=0.35, sE2=0.01, eps=0.5)
GAU = dict(sr=0.1, sp=0.1)

def energy_of_J(J):
    return -1.0/(2.0*(J + c)**2)

def king_feq(J):
    J = np.asarray(J, float)
    out = np.exp(-energy_of_J(J)/KIN['sE2']) - np.exp(-energy_of_J(KIN['Jt'])/KIN['sE2'])
    return np.where(J >= KIN['Jt'], 0.0, out)

def df0(Q, J, dftype):
    Q = np.asarray(Q, float); J = np.asarray(J, float)
    if dftype == 'gauss':
        return np.exp(-np.sin(0.5*Q)**2/GAU['sp']**2)*np.exp(-J**2/GAU['sr']**2)*J**2
    if dftype == 'bimodal':
        p = BIM
        return ((1.0 + p['b1']*np.cos(Q) + p['b2']*np.cos(2*Q))
                *(np.exp(-(J-p['Ja'])**2/p['sa']**2)
                  + p['wb']*np.exp(-(J-p['Jb'])**2/p['sb']**2)))
    if dftype == 'spiral':
        p = SPI
        return (np.exp(-(J-p['J0'])**2/p['s0']**2)
                *np.exp(-np.sin(0.5*(Q - p['beta']*J))**2/p['sq']**2))
    if dftype == 'king':
        return king_feq(J)*(1.0 + KIN['eps']*np.cos(Q))
    raise ValueError(dftype)

def Jrange(dftype):
    return {'gauss':(1e-4*GAU['sr'], 6*GAU['sr']),
            'bimodal':(1e-4, 0.45), 'spiral':(1e-4, 0.45),
            'king':(1e-4, 1.03*KIN['Jt'])}[dftype]

def rp_to_QJ(r, p, L=L0):
    """(r,p_r) -> (Q3,J3), el mismo mapeo que analysish.f90."""
    E  = -1.0/(1.0+np.sqrt(1.0+r**2)) + 0.5*L**2/r**2 + 0.5*p**2
    rt = np.sqrt(1.0 + 2.0*E*(2.0 + 2.0*E + L**2))
    er1 = np.sqrt((1.0 + E*(2.0+L**2) - rt)/(2.0*E**2))
    er2 = np.sqrt((1.0 + E*(2.0+L**2) + rt)/(2.0*E**2))
    s1 = 1.0+np.sqrt(1.0+er1**2); s2 = 1.0+np.sqrt(1.0+er2**2)
    ss = 1.0+np.sqrt(1.0+r**2)
    arg = (s1+s2-2.0*ss)/(s2-s1)
    J = 1.0/np.sqrt(-2.0*E) - c
    cl = np.sign(arg)*np.minimum(np.abs(arg), 1.0)
    eta = np.where(p >= 0.0, np.arccos(cl), np.arccos(-cl)+np.pi)
    ecc = np.sqrt((-2.0*E)**3)*np.sqrt(-L**2-2.0*E-2.0-0.5/E)/(-2.0*E)
    return np.mod(eta - ecc*np.sin(eta), 2*np.pi), J

def omega(J, L=L0):
    return 1.0/(J + 0.5*(L+np.sqrt(L**2+4.0)))**3
