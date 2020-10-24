import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time 
start = time.time()

X_ini = float(input('X initial : '))
Y_ini = float(input('Y initial : '))
Z_ini = float(input('Z initial : '))

b=5.0
X = []
Y = []
Z = []
X.append(X_ini)
Y.append(Y_ini)
Z.append(Z_ini)

#defining the integral params
tstart=0
tend=50000.0
dt = 0.1
t=np.arange(tstart, tend, dt)
N = int((tend - tstart)/dt)


def func_dot(arr):
    der_arr = np.array([arr[1], -arr[0]-arr[2]*arr[1], b*((arr[1]*arr[1])-1)])
    return(der_arr)

'''
def euler(det_func, arr):
    der_arr = det_func(arr)
    delta_val_arr = der_arr*dt
    arr_np1 = np.add(np.asarray(arr),np.asarray(delta_val_arr))
    print(np.shape(arr_np1))
    return(arr_np1)
'''
def rk4(der_func, arr):
    np.asarray(arr)
    k1 = np.array(der_func(arr))*dt
    k2_arg = np.add(arr, k1/2)
    k2 = der_func(k2_arg)*dt
    k3_arg = np.add(arr, k2/2)
    k3 = der_func(k3_arg)*dt
    k4_arg = np.add(arr, k3)
    k4 = der_func(k4_arg)*dt
    delta_arr_1 = np.add(k1, 2*k2)
    delta_arr_2 = np.add(2*k3, k4)
    delta_val_arr = np.add(delta_arr_1, delta_arr_2)*(1/6)
    arr_np1 = np.add(delta_val_arr, arr)
    return(arr_np1)


arr_n = [X_ini, Y_ini, Z_ini]
X_cros = []
Y_cros = []
lyapunov_exp_x = []
lyapunov_exp_y = []
lyapunov_exp_z = []
array_0 = rk4(func_dot, arr_n)
x_0 = array_0[0]
y_0 = array_0[1]
z_0 = array_0[2]

#loop for finding traj
for i in range(N):
    array_np1 = rk4(func_dot, arr_n)
    x_1 = array_np1[0]
    y_1 = array_np1[1]
    z_1 = array_np1[2]
    X.append(array_np1[0])
    Y.append(array_np1[1])
    Z.append(array_np1[2])
    lyapunov_exp_x.append((1/float(i+1))*np.log(np.abs(x_1/x_0)))
    lyapunov_exp_y.append((1/float(i+1))*np.log(np.abs(y_1/y_0)))
    lyapunov_exp_z.append((1/float(i+1))*np.log(np.abs(z_1/z_0)))
    arr_n = array_np1
    if array_np1[2] < 0.01 and array_np1[2] > -0.01:
        X_cros.append(array_np1[0])
        Y_cros.append(array_np1[1])
    if i==int(N/2):
        print('half')
print(X[N-1], Y[N-1], Z[N-1])
'''
#print(time.time()-start)
plt.plot(t,X[1:])
plt.show()
plt.clf()
plt.plot(t,Y[1:])
plt.show()
plt.clf()
plt.plot(t,Z[1:])
plt.show()
plt.clf()
'''
fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.plot(X,Y,Z)
plt.xlabel("X")
plt.ylabel("Y")
#plt.zlabel("Z")
plt.show()

plt.plot(X_cros, Y_cros, '.')
plt.show()
plt.clf()
'''
plt.plot(t, lyapunov_exp_x)
plt.show()
plt.clf()
plt.plot(t, lyapunov_exp_y)
plt.show()
plt.clf()
plt.plot(t, lyapunov_exp_z)
plt.show()
plt.clf()
'''
