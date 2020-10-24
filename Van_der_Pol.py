import math as mt
import numpy as np
import matplotlib.pyplot as plt

# Van der Pol Oscillator
# Hah, lol
t_start=0.
t_end=20.
N=1000
v_start=2.
x_start=2.
#t_start = int(input('Initial Time: '))
#t_end = int(input('Final Time: '))
#N = int(input('Precision: '))
dt = float(t_end-t_start)/float(N)
#x_start = int(input('Initial Position: '))
#v_start = int(input('Initial Velocity: '))
mu = float(input('Strength of Damping: '))

# Runge-Kutta-Nystrom Method for General Second Order ODE
# RKN5 ftw, y'all!

T=np.arange(t_start, t_end+dt, dt)

X = np.zeros(N+1)
X[0] = x_start
V = np.zeros(N+1)
V[0] = v_start

k1=np.zeros(N)
k2=np.zeros(N)
k3=np.zeros(N)
k4=np.zeros(N)
l1=np.zeros(N)
l2=np.zeros(N)
l3=np.zeros(N)
l4=np.zeros(N)

def f(t,x,v):
	return v

def g(t,x,v):
	return (-x + mu*(1-x**2)*v)

for i in range(0,N):
	k1[i]= dt*f(T[i], X[i], V[i])
	l1[i]= dt*g(T[i], X[i], V[i])
	k2[i]=dt*f(T[i]+0.5*dt, X[i]+0.5*k1[i], V[i]+0.5*l1[i])
	l2[i]=dt*g(T[i]+0.5*dt, X[i]+0.5*k1[i], V[i]+0.5*l1[i])
	k3[i]=dt*f(T[i]+0.5*dt, X[i]+0.5*k2[i], V[i]+0.5*l2[i])
	l3[i]=dt*g(T[i]+0.5*dt, X[i]+0.5*k2[i], V[i]+0.5*l2[i])
	k4[i]=dt*f(T[i]+dt, X[i]+k3[i], V[i]+l3[i])
	l4[i]=dt*g(T[i]+dt, X[i]+k3[i], V[i]+l3[i])
	X[i+1]=X[i]+(1.0/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
	V[i+1]=V[i]+(1.0/6.0)*(l1[i]+2*l2[i]+2*l3[i]+l4[i])

plt.subplot(3,1,1)
plt.plot(T,X, 'r', lw=2.0)
plt.xlabel('Time')
plt.ylabel('Position')
plt.grid()
plt.legend()

plt.subplot(3,1,2)
plt.plot(T,V, 'k', lw=2.0)
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.grid()

plt.subplot(3,1,3)
plt.plot(X,V, 'b', lw=2.0)
plt.xlabel('Position')
plt.ylabel('Velocity')
plt.grid()

plt.show()

