from numpy import *
from matplotlib.pyplot import *
def a(t):
	return 1.0

def b(t, c, u0):
	return c + a(t)*(c*t+u0)	

def f_linear(un, tn, u0, c=1.23):
	return -a(tn)*un + b(tn, c, u0)
def f_exp(un, tn, u0):
	return -un + 1
def analytical_exp(t):
	return 1-exp(-t)

def forwardEuler(f,tn,un, dt):
	return un + dt*f(un, tn, un)

def solve(f, u0, t_start, t_stop, N):

	dt = (t_stop - t_start)/(N-1)
	print dt
	#need to find u1. Does this by Forward Euler
	u1 = forwardEuler(f, t_start, u0, dt)
	u = zeros(N)
	t = linspace(t_start, t_stop, N)
	u[0] = u0
	u[1] = u1
	for i in xrange(1,N-1):
		u[i+1] =  2*dt*f(u[i], t[i], u0) + u[i-1]
		
	return t, u


if __name__=='__main__':
	u0=0.1
	t_start = 0
	t_stop = 4.0
	N = 42
	t_l,u_l = solve(f_linear, u0, t_start, t_stop, N)
	plot(t_l,u_l)
	figure()
	t_e,u_e = solve(f_exp, 0, t_start, t_stop, N)
	plot(t_e, u_e)
	hold('on')
	plot(t_e, analytical_exp(t_e))
	legend(["leapfrog", "analytical"])
	show()
	
	
