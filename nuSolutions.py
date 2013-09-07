import numpy as np
import ROOT as r
import math
try: from scipy.optimize import leastsq
except: leastsq=None

def UnitCircle():
    '''Unit circle in extended representation.'''
    return np.diag([1,1,-1])


def cofactor(A,(i,j)):
    '''Cofactor[i,j] of 3x3 matrix A.'''
    a = A[not i:2 if i==2 else None:2 if i==1 else 1,
          not j:2 if j==2 else None:2 if j==1 else 1]
    return (-1)**(i+j) * (a[0,0]*a[1,1] - a[1,0]*a[0,1])


def R(axis, angle):
    '''Rotation matrix about x,y, or z axis (0,1, or 2, respectively.)'''
    c,s = math.cos(angle),math.sin(angle)
    R = c * np.eye(3)
    for i in [-1,0,1] : R[ (axis-i)%3, (axis+i)%3 ] = i*s + (1 - i*i)
    return R


def Derivative():
    '''Matrix to differentiate [cos(t),sin(t),1] and R(2,_).'''
    return R(2, math.pi / 2).dot(np.diag([1, 1, 0]))


def multisqrt(y):
    '''Valid real solutions to y=x*x.'''
    return [] if y<0 else [0] if y==0 else (lambda r: [-r,r])(math.sqrt(y))


def factor_degenerate(G, zero=0):
    '''Linear factors of a degenerate quadratic polynomial.'''
    if G[0,0] == 0 == G[1,1]:
        return [[G[0,1], 0, G[1,2]],
                [0, G[0,1], G[0,2] - G[1,2]]]

    swapXY = abs(G[0,0]) > abs(G[1,1])
    Q = G[(1,0,2),][:,(1,0,2)] if swapXY else G
    Q /= Q[1,1]
    q22 = cofactor(Q,(2,2))

    if -q22 <= zero:
        lines = [[Q[0,1], Q[1,1], Q[1,2]+s] for s in multisqrt(-cofactor(Q,(0,0)))]
    else:
        x0,y0 = [cofactor(Q,(i,2)) / q22 for i in [0,1]]
        lines = [[m, Q[1,1], -Q[1,1]*y0 - m*x0] for m in [Q[0,1] + s for s in multisqrt(-q22)]]

    return [[L[swapXY],L[not swapXY],L[2]] for L in lines]


def intersections_ellipse_line(ellipse=None, line=None, zero=1e-12):
    '''Points of intersection between ellipse and line.'''
    _,V = np.linalg.eig(np.cross(line,ellipse).T)
    sols =  sorted([(v.real / v[2].real,
                     np.dot(line,v.real)**2 +
                     np.dot(v.real,ellipse).dot(v.real)**2)
                    for v in V.T ],
                   key=lambda (s,k): k)[:2]
    return [s for s,k in sols if k<zero]


def intersections_ellipses(A, B, returnLines=False):
    '''Points of intersection between two ellipses.'''
    if abs(np.linalg.det(B)) > abs(np.linalg.det(A)): A,B = B,A
    eigval = next(e.real for e in np.linalg.eigvals(np.linalg.inv(A).dot(B)) if not e.imag)
    lines = factor_degenerate(B - eigval*A)
    points = sum([intersections_ellipse_line(A,L) for L in lines],[])
    return (points,lines) if returnLines else points


class nuSolutionSet(object):
    '''Collection of definitions for neutrino analytic solution in t->b,mu,nu.'''
    
    def __init__(self, b, mu, # Lorentz Vectors
                 Wm2=80.385**2, Tm2=172.5**2 #GeV**2
                 ):
        self.mu = mu
        self.b = b
        self.D2 = 0.5 * (Tm2 - Wm2 - b.M2())
        self.c = r.Math.VectorUtil.CosTheta(b,mu)
        self.s = math.sqrt(1-self.c**2)
        self.x0 = -0.5 * Wm2 / mu.P()
        self.y0 = - ( self.x0*self.c + self.D2 / b.P() ) / self.s

        self.m = (math.sqrt(self.b.E()**2)/self.b.P() - self.c) / self.s
        self.m_alt = (-math.sqrt(self.b.E()**2)/self.b.P() - self.c) / self.s
        self.y1 = -self.x0 / self.m
        self.x1 = self.x0 + (self.y1 - self.y0) / self.m
        self.Z = math.sqrt( max( 0,
                                 self.y1 * (self.y1 - 2*self.y0) - self.x0**2 - Wm2 ) )
        self.Wm2 = Wm2

    @property
    def K(self):
        '''Extended rotation from F' to F coordinates.'''
        return np.array([[self.c, -self.s, 0, 0],
                         [self.s,  self.c, 0, 0],
                         [     0,       0, 1, 0],
                         [     0,       0, 0, 1]])

    @property
    def A_mu(self):
        '''F coordinates constraint on W momentum: paraboloid.'''
        return np.array([[0,       0, 0,                self.x0],
                         [0,       1, 0,                      0],
                         [0,       0, 1,                      0],
                         [self.x0, 0, 0, self.Wm2 - self.x0**2]])

    @property
    def A_b(self):
        '''F coordinates constraint on W momentum: ellipsoid.'''
        K, D2, Wm2, e2 = self.K, self.D2, self.Wm2, self.b.E()**2
        return K.dot(np.array([[ self.b.M2(),   0,  0,  -D2*self.b.P()],
                               [           0,   e2, 0,               0],
                               [           0,   0, e2,               0],
                               [-D2*self.b.P(), 0, 0, e2 * Wm2 - D2**2]
                               ]
                              ) ).dot(K.T)

    @property
    def R_T(self):
        '''Rotation from F coordinates to laboratory coordinates.'''
        R_z = R(2, -self.mu.phi() )
        R_y = R(1,  0.5*math.pi - self.mu.theta() )
        R_x = next( R(0,-math.atan2(z,y)) for x,y,z in (R_y.dot( R_z.dot( [self.b.x(),
                                                                           self.b.y(),
                                                                           self.b.z()]) ),))
        return R_z.T.dot(R_y.T.dot(R_x.T))

    @property
    def H_tilde(self):
        '''Transformation of t=[c,s,1] to nu momentum: F coordinates.'''
        Z, m, x1, y1, p = self.Z, self.m, self.x1, self.y1, self.mu.P()
        return np.array([[ Z/m,  0,  x1 - p ],
                         [ Z,    0,      y1 ],
                         [ 0,    Z,       0 ]])
    @property
    def H(self):
        '''Transformation of t=[c,s,1] to nu momentum: lab coordinates.'''
        return self.R_T.dot(self.H_tilde)

    @property
    def H_perp(self):
        '''Transformation of t=[c,s,1] to nu transverse momentum: lab coordinates.'''
        return np.vstack([self.H[:2],[0,0,1]])
    
    @property
    def N(self):
        '''Solution ellipse of nu transverse momentum: lab coordinates.'''
        E = self.H_perp
        return np.linalg.inv(E.T).dot(UnitCircle()).dot(np.linalg.inv(E))


class singleNeutrinoSolution(object):
    '''Calculation of most compatible neutrino momentum for tt -> lepton+jets.'''
    def __init__(self, b, mu, # Lorentz Vectors
                 (metX,metY), #GeV
                 sigma2, #[GeV**2]
                 Wm2=80.385**2, Tm2=172.5**2 #GeV**2
                 ):
        self.solutionSet = nuSolutionSet(b, mu, Wm2, Tm2)
        S2 = np.vstack( [np.vstack( [ np.linalg.inv(sigma2), [0,0]] ).T, [0,0,0]] )
        V0 = np.outer( [metX, metY, 0 ], [0,0,1])
        deltaNu = V0 - self.solutionSet.H
        
        self.X = np.dot(deltaNu.T, S2).dot(deltaNu)
        M = next( XD + XD.T for XD in (self.X.dot(Derivative()),))
        
        self.solutions = sorted(intersections_ellipses(M,UnitCircle()), key=self.calcX2)
        
    def calcX2(self,t): 
        return np.dot(t,self.X).dot(t)

    @property
    def chi2(self):
        return self.calcX2(self.solutions[0])

    @property
    def nu(self):
        '''Solution for neutrino momentum.'''
        return self.solutionSet.H.dot(self.solutions[0])


class doubleNeutrinoSolutions(object):
    '''Calculation of solution pairs of neutrino momentum, tt -> leptons.'''
    def __init__(self, (b,b_), (mu,mu_), # Lorentz Vectors
                 (metX,metY), #(GeV,GeV)
                 Wm2=80.385**2, Tm2=172.5**2 #GeV**2
                 ):
        self.solutionSets = [nuSolutionSet(B, M, Wm2, Tm2) for B,M in zip((b,b_),(mu,mu_))]

        V0 = np.outer( [metX, metY, 0 ], [0,0,1])
        self.S = V0 - UnitCircle()

        N,N_ = [ss.N for ss in self.solutionSets]
        n_ = self.S.T.dot( N_ ).dot(self.S)

        v = intersections_ellipses(N, n_)
        v_ = [self.S.dot(sol) for sol in v]

        if not v and leastsq:
            es = [ss.H_perp for ss in self.solutionSets]
            met = np.array([metX,metY,1])
            def nus(ts) : return tuple(e.dot([math.cos(t),math.sin(t),1]) for e,t in zip(es,ts))
            ts,_ = leastsq( lambda params : sum( nus(params), -met)[:2], [0,0], ftol=5e-5, epsfcn=0.01 )
            v,v_ = [[i] for i in nus(ts)]

        self.perp = v
        self.perp_ = v_
        self.n_ = n_

    @property
    def nunu_s(self):
        '''Solution pairs for neutrino momenta.'''
        K,K_ = [ss.H.dot(np.linalg.inv(ss.H_perp)) for ss in self.solutionSets]
        return [(K.dot(s), K_.dot(s_)) for s,s_ in zip(self.perp,self.perp_)]
