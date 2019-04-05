#coding=utf-8
'''
计物大作业第七题
求解定态氢原子波函数
'''

import math as ms
import matplotlib.pyplot as plt
import numpy as np # 这两个库仅用作画图用

def LaguerreL(n,a,x):
    '''
    To calculate generalized Laguerre polinomial 
    return L^a_n(x)
    '''
    Li = 1
    Lip = 1-x+a 
    if n == 1:
        return Lip
    elif n==0:
        return Li
    else :
        for i in range(1,n):
            Li,Lip = Lip , 1/(i+1)*((2*i+1+a-x)*Lip-(i+a)*Li)
    return Lip

def Harmonic(l,m,theta,phi):
    '''
    To calculate Spherical Harmonic function
    return Y_lm(\theta,\phi)
    '''

    x = ms.cos(theta)
    P = GeneralizedLegrend(l,m,x)
    # 直接得到Y
    Y = (-1)**m*ms.sqrt((2*l+1)*ms.factorial(l-abs(m))/(4*ms.pi*ms.factorial(l+abs(m))))*P*complex(ms.cos(phi),ms.sin(phi))
    return Y

def GeneralizedLegrend(l,m,x):
    '''
    求解连带勒让德函数
    基本思路是利用递推公式达到二维平面上(l,m)这个点
    先达到(l,0)，求解勒让德函数
    返回连带勒让德P^l_m(x)
    '''
    Pi = 1
    Pip = x
    sq = ms.sqrt(1-x**2) #减少开方运算
    for i in range(1,l):
        Pip , Pi = ((2*i+1)*x*Pip-i*Pi)/(i+1) , Pip
    # 利用(l,0)和(l-1,0)求出(l,1)
    Pii = l*(x*Pip-Pi)/sq
    # 继续迭代得到(l,m)
    for i in range(1,m):
        Pii , Pip = (-2*i*x*Pii/sq-(l+i)*(l-i+1)*Pip) , Pii
    return Pii

def val(n,l,m,x,y):
    '''
    由于在画图过程中，最外部的归一化系数不需要重复计算，那么直接计算为归一化的函数会比较节约时间
    '''
    rho = 2*ms.sqrt(x**2+y**2)/n
    return ms.exp(-rho)*rho**(2*l)*LaguerreL(n-l-1,2*l+1,rho)**2*(Harmonic(l,m,ms.pi/2,0).real)**2

def draw(n,l,m,r):
    '''
    画图函数,输入波函数Psi_nlm,以及最大绘图范围r
    '''
    max = 10**3
    x = y = np.linspace(-r,r,max)
    z = np.array([val(n,l,m,i,j) for i in x for j in y])
    Z = z.reshape(max,max)
    im = plt.imshow(Z\
        ,extent=[-Z.shape[1]/2., Z.shape[1]/2., -Z.shape[0]/2., Z.shape[0]/2. ])
    plt.colorbar(im)
    plt.show()

    return

if __name__ == '__main__':
    # print(LaguerreL(30,40,100))
    # print(Harmonic(2,1,ms.pi/3,-ms.pi/3))
    # print(GeneralizedLegrend(2,1,0.5))
    # draw(2,1,0,5)
    