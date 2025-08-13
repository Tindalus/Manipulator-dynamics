import numpy as np
from sympy.physics.vector import dynamicsymbols
import sympy as sp
from sympy import *

def LDynamics(DH: list, mass: list, lengths: list, orients: list, Inertions: list, CreateTxt: bool = False):
        '''
        This function returns the SymPy simbol matrixes of n actuators linked chain
        Args:
             DH (listlike (n,4)): DH table of your system's reference frame parameters
             mass (list (n)): list of all link and tool masses
             lengths (list (n)): list of all link and tool lengths
             orients (list (n)): list of each link's orientation where 0 = z, 1 = y, 2 = x, so for eg. [0,2,1] for 3 link chain.
             Inertions (list (n)): list of each link inertias along main axis
             CreateTxt (bool (optional)): check for create 3 .txt files for each output matrix. Deafault False
        Returns:
             \n
             List (n,n): Inertia Matrix M \n
             List (n,1): Coriolis matrix multiplied on qdot vector\n
             List (n,1): Gravity vector
        '''

        # symbols assignment
        t = sp.symbols('t')
        n = len(mass) #system's degree
        x,y,z,q,qdot,qdotdot,Tau = [0]*n,[0]*n,[0]*n,[0]*n,[0]*n,[0]*n,[0]*n
        
        for i in range(len(x)):
                x[i],y[i],z[i],q[i],qdot[i],qdotdot[i],Tau[i] = sp.symbols(f'x{i},y{i},z{i},q{i},qdot{i},qdotdot{i},Tau{i}')
                q[i] = dynamicsymbols(f'q{i}')
                qdot[i] = dynamicsymbols(f'qdot{i}')

        #physical parameters
        g = 9.81
        m = mass
        l = lengths

        def Ti(DH):
                #for transformation matrixes
                def Ai(i):
                        A = [[cos(DH[i][3]),-1*cos(DH[i][1])*sin(DH[i][3]),sin(DH[i][1])*sin(DH[i][3]),DH[i][0]*cos(DH[i][3])],
                        [sin(DH[i][3]),cos(DH[i][1])*cos(DH[i][3]),-1*sin(DH[i][1])*cos(DH[i][3]),DH[i][0]*sin(DH[i][3])],
                        [0,sin(DH[i][1]),cos(DH[i][1]),DH[i][2]],
                        [0,0,0,1]]
                        return A
                
                A_ = np.eye(4)
                T = []
                for i in range(len(DH)):
                        A_ = np.matmul(A_,Ai(i))
                        T.append(A_)
                return T
        


        # #tranlation of x,y,z through q
        T = Ti(DH)
        for i in range(len(q)):
                x[i],y[i],z[i] = T[i][0][3],T[i][1][3],T[i][2][3]

        #calculating velosity vectors for kinetic energy
        v = []
        for j in range(len(q)):
                v.append([sp.diff(x[j],t),sp.diff(y[j],t),sp.diff(z[j],t)])

        #fixing Derivatives in the output v
        V = []
        for i in range(len(q)):
                V.append([0 for i in range(3)])
                for k in range(3):
                        V[i][k] = v[i][k].subs([(sp.Derivative(q[j], t), qdot[j]) for j in range(len(q))])
        v = V

        #makin angular velocities
        Ii = Inertions
        w = []
        W1,W2,W3 = 0,0,0
        for i in range(len(q)):
                if orients[i] == 0:
                        W1 += qdot[i]
                        w.append(W1)
                elif orients[i] == 1:
                        W2 += qdot[i]
                        w.append(W2)
                elif orients[i] == 2:
                        W3 += qdot[i]
                        w.append(W3)

        #Forming kinetic energy for lagrangian
        K1 = 0
        for i in range(len(v)):
                vtv = np.matmul(np.transpose(v[i]),v[i])
                K1 += m[i]*vtv
        K1 *= 1/2
        K2 = 0
        for i in range(len(q)):
                K2 += w[i]*w[i]*Ii[i]
        K2 *= 1/2
        K = K1 + K2 #kinetic energy
        
        #calculating potential energy
        #taking -x as a height for potential energy
        P = 0
        for i in range(len(q)):
                if n <= 2:
                        P += m[i]*g*(1*y[i]) #choose y as gravity axis by default when 2d 
                else:
                        P += m[i]*g*(1*z[i])

        #Lagrangian
        L = K - P

        #finding Lagrange equation left side
        LagrangeLeft = [] #list containing left sides of the Lagrange equation = 0
        for i in range(len(q)):
                LagrangeLeft.append(sp.diff(sp.diff(L,qdot[i]),t) - sp.diff(L,q[i])  ) # can put here  - Tau[i]
        #fixing Derivatives in the output LagrangeLeft
        V = []
        sub = []
        for j in range(len(q)):
                sub.append(((sp.Derivative(q[j], t), qdot[j])))
                sub.append((sp.Derivative(qdot[j], t), qdotdot[j]))
        for i in range(len(q)):
                V.append(0)
                V[i] = LagrangeLeft[i].subs(sub)
        LagrangeLeft = V



        #find all quadratic coefficients in K and form matrix out of them
        #It is an M(q) inertia matrix in fact
        M = []
        for i in range(len(q)):#rows
                M.append([0 for i in range(len(q))])
                for j in range(len(q)):#collums
                        M[i][j] = (sp.expand(K)).coeff(qdot[i]*qdot[j],1)*2


        #C(q,qdot)*qdot matrix
        Mdot = []
        for i in range(len(q)):
                Mdot.append([0]*n)
                for j in range(len(q)):
                        Mdot[i][j] = sp.diff(M[i][j],t)
        C_qdot1 = np.matmul(Mdot,qdot)
        C_qdot2 = []
        for i in range(len(q)):
                C_qdot2.append(sp.diff(np.matmul(np.transpose(qdot),np.matmul(M,qdot)),q[i]))
        C_qdot = [C_qdot1[j]-0.5*C_qdot2[j] for j in range(len(q))]
        #fixing Derivatives in the output C_qdot
        V = []
        for i in range(len(q)):
                V.append(0)
                V[i] = C_qdot[i].subs([(sp.Derivative(q[j], t), qdot[j]) for j in range(len(q))])
        C_qdot = V


        #g(q) gravity term
        G = []
        for i in range(len(q)):
                G.append(sp.diff(P,q[i]))
        #fixing Derivatives in the output G
        V = []
        for i in range(len(q)):
                V.append(0)
                V[i] = G[i].subs([(sp.Derivative(q[j], t), qdot[j]) for j in range(len(q))])
        G = V


        x,xdot = [0]*n,[0]*n
        for i in range(len(x)):
                x[i],xdot[i] = sp.symbols(f'x{i},xdot{i}')
        #fixing (t) in the outputs
        sub = []
        for j in range(len(q)):
                sub.append(((q[j], x[j])))
                sub.append((qdot[j], xdot[j]))
        V = []
        for i in range(len(x)):
                V.append([0]*n)
                for j in range(len(x)):
                        V[i][j] = M[i][j].subs(sub)
        M = V
        V = []
        for i in range(len(x)):
                V.append(0)
                V[i] = C_qdot[i].subs(sub)
        C_qdot = V
        V = []
        for i in range(len(x)):
                V.append(0)
                V[i] = G[i].subs(sub)
        G = V

        # print("M(q)=")
        # # for i in range(len(M)):
        # #         print(M[i])
        # print(M)
        # print("C(q,qdot)*qdot=")
        # print(C_qdot)
        # print("G(q)=")
        # print(G)

        #saving output as .txt files
        if CreateTxt == 1:
                f = open("M.txt", "w")
                f.write(str(M))
                f.close()
                f = open("C.txt", "w")
                f.write(str(C_qdot))
                f.close()
                f = open("G.txt", "w")
                f.write(str(G))
                f.close()
        return M,C_qdot,G
