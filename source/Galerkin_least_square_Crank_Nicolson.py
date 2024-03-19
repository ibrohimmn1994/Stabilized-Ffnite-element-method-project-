from dolfin import*
import numpy as np
import scipy.linalg as la
import ufl



T = 2
num_samples = 40
dt = 0.05


mesh = RectangleMesh(Point(-2,-2.5),Point(2,1.5),31,31)
V = FunctionSpace(mesh, "CG", 1)
VV = VectorFunctionSpace(mesh,"CG",1)

def boundary(x, on_boundary):
    return on_boundary


bc = DirichletBC(V,Constant(0.0),boundary)

ID = Expression('pow(pow(x[0],2)+pow(x[1],2),0.5) <=1?14/4*pi:pi/4', degree=2,pi=np.pi)


u0 = interpolate(ID, V)

u = TrialFunction(V)
v = TestFunction(V)
#B = Expression(('sin(u)','cos(u)'),degree=2, u=u)
#B0 = Expression(('sin(u)','cos(u)'), degree=2,u=u0)

def B(u):
   # return project(Expression(('-cos(u)','sin(u)'),degree=2,u=u),VV )
    return    as_vector(( cos(u), -sin(u) ))  
  # return np.array([sym(cos(u)),sym(sin(u))])
    

#def f(u):
 #   return as_vector((sin(u),cos(u)))

u_k = interpolate(Constant(0.0),V)
eps=1.0
tol = 1.0E-5
maxiter=25

h = 4/31
gamma = 0.5*h

a = u*v*dx + 0.5*dt*dot(B(u_k),grad(u) ) *v*dx\
        +gamma*u*dot(B(u0),grad(v))*dx \
        +0.5*gamma*dt*dot(B(u_k),grad(u))*dot(B(u0),grad(v))*dx

L = u0*v*dx - 0.5*dt*dot(B(u0),grad(u0) )*v*dx \
        +gamma*u0*dot(B(u0),grad(v))*dx\
        -0.5*gamma*dt*dot(B(u0),grad(u0))*dot(B(u0),grad(v))*dx

u = Function(V)
 


out_file = File("VTK/Results.pvd", "compressed")


#u.assign(u0)
t=0
t_save =0.0
out_file << (u0,t)

while t<=T:
    t += dt
    t_save += dt

    itera = 0
    eps =np.Inf

    while eps>tol and itera <maxiter:
        itera += 1
        
        solve(a==L,u,bc)
        diff = np.array(u.vector()) - np.array(u_k.vector())
        
        eps = np.linalg.norm(diff, ord=np.Inf)
        u_k.assign(u)
   
    u0.assign(u)
    
    if t_save >T/num_samples or t>=T-dt:
        print('time = ',t)
        out_file << (u,t)
        t_save =0
























