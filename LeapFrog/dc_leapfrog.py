from numpy import *
from matplotlib.pyplot import *

class Model:
    
    def __init__(self,a=lambda t:t, b=None,c=1, I= 0.1, analytic = None):
        """
        The default model has a linear solution. To avoid linear solution you have to provide b not None! a(t) and b(t) and analytic(t) should be functions, c and I a constant. I is the inital condition u(0) = I.
        """
        
        self.a = a
        self.c = c
        self.I = I
        if b==None:
            self.b = lambda t: c + a(t)*(c*t+I)   #now we know the sol is linear
            self.analytic = self.analyticalLinearSol #set analytical sol -->linear
                                                 
        else:
            self.b = b
            self.analytic = analytic # b i given, so we need to know the analytic solution(cannot assume linear solution anymore, it might be linear but we dont know)

  
    def f(self,un,tn):
        """
        The general equation u'(t) = f(u,t), in our model u'(t) = -a(t)u(t) + b(t)
        """
        return -self.a(tn)*un + self.b(tn)
    
    def analyticalLinearSol(self, t):
        """
        Only for equations with linear solutions
        """
        return self.c*t + self.I


class Solver:

    def __init__(self, model, t_start=0, t_stop=4.0, N=100):
        """
        f should be a function that takes to arguments f(u,t), u0 the initial condition, t_start/t_stop defines the interval, N is number of mesh points.
        """
        self.model = model
        self.f = model.f
        self.u0 = model.I
        self.t_start = t_start
        self.t_stop = t_stop
        self.N = N
        self.dt = (t_stop - t_start)/(N-1)
        self.u1 = self.forwardEuler(self.u0, t_start)
        
        self.t = linspace(t_start, t_stop, N)
       

    def forwardEuler(self,un, tn):
        """
        Used to find u1
        """
        return un + self.dt*self.f(un, tn)

    def solve_LF(self):
        """
        Solves the differential eq with the Leapfrog scheme
        """
        self.u = zeros(self.N)
        self.u[0] = self.u0
        self.u[1] = self.u1
        u = self.u
        f= self.f
        dt = self.dt
        t = self.t
        N = self.N
        for n in xrange(1,N-1):
            u[n+1] = 2*dt*f(u[n],t[n]) + u[n-1]
        #return t,u
    
    def plot_num_analy_sol(self):
        """
        Plots the numerical solution, and the analytical if the solution is linear, or the analytical solution is given.
        """
        figure()
        plot(self.t,self.u)
        xlabel("time (s)")
        ylabel("u(t)")
        title("Solution to the differential equation u'(t) = -a(t)u(t) +b(t)")
        if self.model.analytic != None:
            hold('on')
            plot(self.t, self.model.analytic(self.t))
            legend(["Numerical solution", "Analytical solution"])
        
        else:
            print "Analytical solution not given"
        show()

    def find_error(self):
        E = None
        if self.model.analytic != None:
            u = self.u
            t = self.t
            u_e = self.model.analytic(t)
            e = u_e -u
            E = sqrt(self.dt*sum(e**2))
        return E
            
        
    def convergence_rate(self, dt_max, dt_min, N_dt):
        dt_values = linspace(dt_max, dt_min, N_dt)
        
        self.E_array = zeros(N_dt)
        self.r_array = zeros(N_dt-1)
        E_array = self.E_array
        
        i = 0
        for dt in dt_values:
            self.dt = dt
            self.N = int((self.t_stop - self.t_start)/dt) + 1
            self.t = linspace(self.t_start,self.t_stop, self.N)
            self.solve_LF()
            self.E_array[i] = self.find_error()
            if i>0:
                n = i-1
                self.r_array[n] = log(E_array[i-1]/E_array[i])/log(dt_values[i-1]/dt)
            i+=1
        figure()
        plot(dt_values[:-1], self.r_array)
        title("Convergence rates for different dt")
        xlabel("dt")
        ylabel("Rate")
        show()

  

    
if __name__ == '__main__':
    linear = Model()
    sol = Solver(linear)
    sol.solve_LF()
    #sol.plot_num_analy_sol()

    def a(t):
        return 1
    def b(t):
        return 1
    u0 = 0
    def exact(t):
        return 1-exp(-t)

    exp_model = Model(a,b,I=0,analytic = exact)
    exp_sol = Solver(exp_model, N = 10000)
    exp_sol.solve_LF()
    exp_sol.plot_num_analy_sol()
    #sol.convergence_rate(1, 0.0001, 20)
    
    
    
    
    
    
