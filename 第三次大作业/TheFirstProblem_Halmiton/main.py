'''
第三次大作业第一题
为了求解电离的Hamiltonian
'''

import numpy as np
from Cholesky import cholesky,solve # 调用自己的Cholesky方法
from copy import deepcopy # 深度拷贝方法
import matplotlib.pyplot as plt

def Build_H(xmax,N):
    '''
    构建哈密顿矩阵
    '''
    Dx = 2*xmax/(N-1)
    D1 = np.diag([1/Dx**2-1/np.sqrt((i)**2+2) for i in np.linspace(-xmax,xmax,N)])
    D1 += np.diag([-1/(2*Dx**2) for i in range(N-1)],k=1)+np.diag([-1/(2*Dx**2) for i in range(N-1)],k=-1)
    return D1

def antipower_method(v,N,xmax,x0=0.5,e=10**(-10),Nmax=10**5):
    '''
    反幂法
    初始化矢量v,格点个数N,解空间xmax,初始位移x0=0.48,精度e=10**(-6),最大迭代次数Nmax=10**5
    '''
    Dx = 2*xmax/(N-1)
    A = [[1/Dx**2-1/np.sqrt((i)**2+2)+x0 for i in np.linspace(-xmax,xmax,N)],[-1/(2*Dx**2) for i in range(N-1)]] # 已经加上位移
    # B = np.diag(A[0])+np.diag(A[1], k=1)+np.diag(A[1], k=-1)
    L = cholesky(A)
    n =0
    # u = v.copy()
    u = v[:] # 深复制
    while n < Nmax:
        # v = solve(L,u.copy())
        v = solve(L, u[:])
        # v = np.linalg.inv(B)@u
        a = 1.0/max(v,key = lambda x:abs(x))
        # v = a*v # 规格化
        v = [a*each for each in v]
        delta = max([v[i]-u[i]for i in range(N)],key = lambda x:abs(x)) # v 已经规格化了
        print(str(n)+" "+str(delta))
        if abs(delta) < e:
            break
        # u =v.copy()
        u = v[:]
        n+=1
    else:
        raise("Two many loops "+str(delta))
    return 0.5, v # 返回的是特征值和本征矢
    # return v , a

def shape(psi,xmax):
    '''
    波函数作图函数
    '''
    N = len(psi)
    x = np.linspace(-xmax,xmax,N)
    y = abs(np.array(psi))**2
    plt.plot(x,y)
    plt.show()


def check(A,v,u):
    A = np.diag(A[0])+np.diag(A[1],k=1)+np.diag(A[1],k=-1)
    e = np.linalg.norm(A@np.asarray(v)-np.asarray(u))
    if abs(e)> 10**(-10):
        print(e)
        # print(A)
        # print("u = "+str(u))
        # print("v = "+str(v))
        return False
    else:
        return True

def groud_state(xmax,N):
    '''
    求解基态波函数
    '''
    # H = Build_H(xmax, N)
    # v = np.array([0.5 for i in range(N)])
    v = [0.5 for i in range(N)]
    # u = antipower_method(np.copy(v), N, xmax)
    u = antipower_method(v, N, xmax)
    return u

def multi(H,u):
    '''
    三对角乘法函数
    '''
    n = len(u)
    b = [0 for _ in range(n)]
    b[0] = H[0][0]*u[0]+H[1][0]*u[1]
    b[n-1] = H[0][n-1]*u[n-1]+H[1][n-2]*u[n-2]
    for i in range(1,n-1):
        b[i] = H[0][i]*u[i]+H[1][i]*u[i+1]+H[1][i-1]*u[i-1]
    return b

def Time_evolution(psi0,xmax,N,dt = .05,T_max = 18*2*np.pi):
    '''
    第一题第二问含时演化问题
    '''
    Dx = 2*xmax/(N-1)
    Nt = int(T_max/dt)
    E0 = np.sqrt(10**(16-16)/3.5094448314)
    H0 = [[1/Dx**2-1/np.sqrt((i)**2+2) for i in np.linspace(-xmax,xmax, N)], [-1/(2*Dx**2) for i in range(N-1)]]  # t=0时哈密顿量
    psi = deepcopy(psi0)
    sequence =[]
    sequence.append(1)
    for t in np.linspace(0,T_max,Nt):
        print("t = "+str(t))
        t1 = t+ .5*dt
        dH = [x*E0*np.sin(t1/(2*18))**2*np.sin(t1) for x in np.linspace(-xmax,xmax,N)] # 计算哈密顿量含时项
        H = [[H0[0][i]+dH[i] for i in range(N)],[each for each in H0[1]]]
        b = multi([[1-.5j*dt*each for each in H[0]],
                   [-.5j*dt*each for each in H[1]]], psi[:])
        L = cholesky([[1+.5*1j*dt*each for each in H[0]],
                      [.5*1j*dt*each for each in H[1]]])
        v = solve(L,deepcopy(b))
        psi = deepcopy(v)
        res = abs(sum([psi[i]*psi0[i].conjugate() for i in range(N)]))**2
        sequence.append(res)
    return sequence,psi

def moment_spec(psi,Dx,dk=0.05):
    N = len(psi)
    k_s = np.linspace(-2.5,2.5,int(5/0.05))
    Pk = []
    for k in k_s:
        Pk.append(Dx/(2*np.pi)**2*abs(sum([np.exp(1j*k*(-N/2+i)*Dx)*psi[i] for i in range(N)]))**2)
    y = np.asarray(Pk)
    # y1 = np.log(y)
    a = plt.subplot(1,2,1)
    plt.plot(k_s,y)
    a.set_title(r"P_k relation")
    # b = plt.subplot(1,2,2)
    # plt.plot(k_s,y1)
    # b.set_title(r"Log(P_k) plot")
    
    plt.show()


def Radiation(psi0, xmax, N, dt=.05, T_max=48*2*np.pi):
    '''
    第一题第四问辐射问题
    '''
    Dx = 2*xmax/(N-1)
    Nt = int(T_max/dt)
    E0 = np.sqrt(2*10**(14-16)/3.5094448314)
    omega = 45.5633525316/300
    H0 = [[1/Dx**2-1/np.sqrt((i)**2+2) for i in np.linspace(-xmax, xmax, N)],
          [-1/(2*Dx**2) for i in range(N-1)]]  # t=0时哈密顿量
    psi = deepcopy(psi0)
    d = []
    a = []
    a.append(0)
    d.append(0)
    for t in np.linspace(0, T_max, Nt):
        print("t = "+str(t))
        t1 = t + .5*dt
        dH = [x*E0*np.sin(omega*t1/(2*48))**2*np.sin(omega*t1)
              for x in np.linspace(-xmax, xmax, N)]  # 计算哈密顿量含时项
        H = [[H0[0][i]+dH[i] for i in range(N)], [each for each in H0[1]]]
        b = multi([[1-.5j*dt*each for each in H[0]],
                   [-.5j*dt*each for each in H[1]]], psi[:])
        L = cholesky([[1+.5*1j*dt*each for each in H[0]],
                      [.5*1j*dt*each for each in H[1]]])
        v = solve(L, deepcopy(b))
        psi = deepcopy(v)
        d1 = sum([abs(psi[i])**2*(Dx*(i-N/2)) for i in range(N)])
        # a1 = sum([abs(psi[i])**2*(E0*np.sin(omega*t/(2*48) )**2*np.sin(omega*t)\
            # -( Dx*(i-N/2) )/( (Dx*(i-N/2))**2 + 2 )**1.5)*Dx for i in range(N)])
        d.append(d1)
        # a.append(a1)
    # 求傅里叶变换
    w = np.linspace(0,15*omega,1000)
    d_omega = [np.log(abs(o**2/np.sqrt(2*np.pi)*sum([np.exp(-1j*o*i*dt)*d[i]*dt for i in range(Nt)]) )**2) for o in w]
    # a_omega=[np.log(abs(1/np.sqrt(2*np.pi)*sum([np.exp(-1j*o*i*dt)*a[i]*dt for i in range(Nt)]) )**2) for o in w]
    x = np.linspace(0,15,1000)
    plt.plot(x,d_omega,label='diplole method')
    # plt.plot(x,a_omega,label='a')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    N = 20000
    xmax=2000
    Dx = 2*xmax/(N-1)
    T_max = 18*2*np.pi
    dt = .05
    Nt = int(T_max/dt)
    # 第一问
    v =np.array([0.5 for i in range(N)])
    lab,psi = antipower_method(v, N, xmax,x0=.50)
    # print(lab)
    # shape(psi,xmax)

    # 第二问
    # psi0 = groud_state(xmax,N)[1]
    # a = np.linalg.norm(np.asarray(psi0))
    # psi0 = [each/a for each in psi0]
    # P,psi = Time_evolution(psi0,xmax,N,dt=dt,T_max=T_max)
    # p0 = sum([psi[i]*psi0[i].conjugate() for i in range(N)])
    # psi_f = [psi[i]-p0*psi0[i] for i in range(N)]
    # moment_spec(psi_f,Dx)
    # t = np.linspace(0,T_max,Nt+1)
    # plt.plot(t,P)
    # plt.show()

    # 第三问
    Radiation(psi,xmax,N,dt=0.05)
