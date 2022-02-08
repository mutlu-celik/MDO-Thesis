# -*- coding: utf-8 -*-
"""
1D Heat conduction
Martin Albers
[1] Strang G. and Fix G. (2008): An analysis of the Finite Element Method,Second Edition, Wellesley-Cambridge Press, Wellesley USA
[2] R. W. Lewis et al. (1996):  The Finite Element Method in Heat Transfer Analysis, John Wiley and Sons, West Sussex England
--- Problem ---
Solving 1-D Heat diffusion in a unit rod
 
"""
import numpy as np
import matplotlib.pyplot as plt

th=np.array([0.008, 0.002])  #thickness
k=np.array([0.5, 14])
rho=np.array([650, 7800])  #density
cp=np.array([1200, 750])  #specific heat
hc=rho*cp #heat capacity
dx=0.001
t_th=np.sum(th) #total thickness of the road
x=np.arange(0, t_th+dx, dx) #spatial discretization
nn=x.size  #nr of nodes
A=1 #area taken as unit

#Boundary conditions
h_mat=np.loadtxt('h.txt')
T_rec_mat=np.loadtxt('T_rec.txt')
time, h=h_mat.T
T_rec=T_rec_mat.T[1][:]
T_rad=250
sb=5.67e-8
e_r=0.8
#time paramaters
dt=0.01
t=time[len(time)-1]
tt=np.arange(0,t+dt,dt)
hi=np.interp(tt,time, h)
T_reci=np.interp(tt,time,T_rec)

e=np.zeros((nn-1,2))
#Creating the element matrix
for i in range(nn-1):
    e[:][i]=np.array([i, i+1])

#Creating space
nel=len(e)
K=np.zeros((nn,nn))
M=np.zeros((nn,nn))
C=np.zeros((nn,nn))
D=np.zeros((nn,nn))
F=np.zeros(nn)
Q=np.zeros(nn)
#c1=np.array([[1, -1], [-1, 1]])
#c2=np.array([[2, 1], [1, 2]])
#Assembling matrices

for j in range(nel):
    nodes=e[:][j]
    if j<(th[0]/dx): #assembling matrices for the first material
       Ke = ((A*k[0])/dx)*np.array([[1, -1], [-1, 1]])
       Me = ((hc[0]*A*dx)/6)*np.array([[2, 1], [1, 2]])
       Fe = ((A*dx)/6)*np.array([[2, 1], [1, 2]])@np.array([Q[int(nodes[0])], Q[int(nodes[1])]])
    for i in range(len(th)):
            if float(j)>=(np.sum(th[0:i+1])/dx) and float(j)<(np.sum(th[0:i+2])/dx):
                Ke = ((A*k[0])/dx)*np.array([[1, -1], [-1, 1]])
                Me = ((hc[0]*A*dx)/6)*np.array([[2, 1], [1, 2]])
                Fe = ((A*dx)/6)*np.array([[2, 1], [1, 2]])@np.array([Q[int(nodes[0])], Q[int(nodes[1])]])   
    K[int(nodes[0]):int(nodes[1])+1,int(nodes[0]):int(nodes[1])+1]=K[int(nodes[0]):int(nodes[1])+1,int(nodes[0]):int(nodes[1])+1]+Ke
    M[int(nodes[0]):int(nodes[1])+1,int(nodes[0]):int(nodes[1])+1]=M[int(nodes[0]):int(nodes[1])+1,int(nodes[0]):int(nodes[1])+1]+Me
    F[int(nodes[0]):int(nodes[1])+1]=F[int(nodes[0]):int(nodes[1])+1]+Fe

#Initial Condition
T_old=np.zeros(nn)+300  
T_new=T_old  

C=(M+0.5*dt*K);             
D=(M-0.5*dt*K);
rr=1
#time stepping
for i in range(int(t/dt+1)):
        K_r=K[len(K)-1][len(K)-1]+hi[i]*A+T_old[nn-1]**3*e_r*sb
        C[len(C)-1][len(C)-1]=M[len(M)-1][len(M)-1]+0.5*dt*K_r
        D[len(D)-1][len(D)-1]=M[len(M)-1][len(M)-1]-0.5*dt*K_r
        F[len(F)-1]=hi[i]*T_reci[i]+T_rad**3*e_r*sb
#        if i%5==0:
        T_new=np.linalg.inv(C)@(D@T_old+dt*F)
        T_old=T_new
        
plt.plot(x,T_new)
plt.show()


