import numpy as np
from numpy.linalg import inv
from LagrangianDynamics import LDynamics
from sympy.physics.vector import dynamicsymbols
from sympy import pi,symbols,Matrix,simplify

#initials
n = 2 #systems degree
Ts = 0.1 #numerical grid step
Tstop = 1 #simulation time
N = int(Tstop/Ts) #number of points

Coords = []
Velosities = []
U = np.zeros(shape=(n,1)) #control input
x,xdot,q,qdot = [0]*n,[0]*n,[0]*n,[0]*n #coordinates and velosities
for i in range(n):
    x[i] = symbols(f'x{i}')
    xdot[i] = symbols(f'xdot{i}')
    q[i] = symbols(f'q{i}')
    q[i] = dynamicsymbols(f'q{i}')
    qdot[i] = dynamicsymbols(f'qdot{i}')


#parameters of model
l1 = 0.25
l2 = 0.22
lengths = [l1, l2]
mass = [0.4,0.35] #mass list
orients = [0,0]
Inertions = [0.00075,0.00075]

#D-H table
#      a        α    d    θ
DH = [[l1,      0,   0, q[0]],                # 0 -> 0
      [l2,      0,   0, q[1]]]                # 0 -> 1

#Dynamics matrixes
M,C,G = LDynamics(DH,mass,lengths,orients,Inertions)
M,C,G = Matrix(M),Matrix(C),Matrix(G)
x,xdot = [0]*n,[0]*n

for j in range(N):
      #Matrix evaluation block
      var_dict = {}
      for i in range(n):
            var_dict[f'x{i}'] = float(x[i])
            var_dict[f'xdot{i}'] = float(xdot[i])
      #evaluate matrixes for current step
      Mevaluated = np.array(M.evalf(subs=var_dict))
      Cevaluated = np.array(C.evalf(subs=var_dict))
      Gevaluated = np.array(G.evalf(subs=var_dict))

      x = Matrix(x)
      x = np.array(x.evalf(subs=var_dict))
      xdot = Matrix(xdot)
      xdot = np.array(xdot.evalf(subs=var_dict))
      
      U = Gevaluated+0.1 #canceled gravity term + some constant signal

      Xdot1 = xdot + Ts*np.matmul(inv(np.array(Mevaluated, dtype=np.float64)), U-Cevaluated-Gevaluated) #new velocity
      Xdotdot = (Xdot1 - xdot)/Ts #acceleration
      xdot = Xdot1 #update
      x1 = x + Ts*xdot #integrate
      x = x1 #update

      #store current step
      Coords.append(x)
      Velosities.append(xdot)


#plotting the result
import matplotlib.pyplot as plt 
fig, axes = plt.subplots(nrows=2, ncols=2)
t = np.arange(0,Tstop,Ts)
axes[0,0].plot(t, np.array(Coords)[:,0], color='red', linewidth=1, markersize=1) 
axes[0,0].set_title('x0(t)')
axes[0,0].grid()
axes[0,0].tick_params(axis='both', which='major', pad=15)
axes[1,0].plot(t, np.array(Coords)[:,1], color='red', linewidth=1, markersize=1) 
axes[1,0].set_title('x1(t)')
axes[1,0].grid()
axes[1,0].tick_params(axis='both', which='major', pad=15)
axes[0,1].plot(t, np.array(Velosities)[:,0], color='red', linewidth=1, markersize=1) 
axes[0,1].set_title('xdot0(t)')
axes[0,1].grid()
axes[0,1].tick_params(axis='both', which='major', pad=15)
axes[1,1].plot(t, np.array(Velosities)[:,1], color='red', linewidth=1, markersize=1) 
axes[1,1].set_title('xdot1(t)')
axes[1,1].grid()
axes[1,1].tick_params(axis='both', which='major', pad=15)

plt.tight_layout()
plt.show()