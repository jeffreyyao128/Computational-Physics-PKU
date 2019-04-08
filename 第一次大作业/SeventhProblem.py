#coding=utf-8
'''
计物大作业第七题
求解定态氢原子波函数
'''

import math as ms # 用于开根、cos计算等
import matplotlib.pyplot as plt 
import numpy as np # 这两个库仅用作画图用
from scipy.special import genlaguerre # 仅用来检验广义拉盖尔函数 注释掉这一行不影响使用

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
    return Y_lm(\ theta,\phi)
    '''

    x = ms.cos(theta)
    P = GeneralizedLegrend(l,m,x)
    # 直接得到Y
    Y = (-1)**m*ms.sqrt((2*l+1)*ms.factorial(l-abs(m))/(4*ms.pi*ms.factorial(l+abs(m))))*P*complex(ms.cos(phi),ms.sin(phi))
    return Y

def test(n,a):
    '''
    检验拉盖尔多项式
    '''
    x = np.linspace(0,10,10)
    y1 = np.array([genlaguerre(n,a)(i)-LaguerreL(n,a,i) for i in range(0,10)])
    plt.plot(x,y1)
    plt.title(r'Error of Laguerre with n= '+str(n)+r' $\alpha = $'+str(a))
    plt.xlabel('x')
    plt.ylabel(r'$\Delta L_n^\alpha(x)$')
    plt.show()

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
    # print(sq)
    # print(Pii)
    # 继续迭代得到(l,m)
    for i in range(1,m):
        Pii , Pip = (-2*i*x*Pii/sq-(l+i)*(l-i+1)*Pip) , Pii
    return Pii

def New_GeneralizedLegrend(l,m,x):
    '''
    新的连带勒让德函数计算方法，减少溢出风险
    '''
    P = 0
    Pi = 1# 并不直接代入P_m^m(x)
    if l-m< 0:
        return 0
    for a in range(m+1,l+1):
        Pi , P = x*ms.sqrt((4*a**2-1.0)/(a**2-m**2))*Pi-ms.sqrt((2*a+1.0)/(2*a-3.0)*((a-1.0)**2-m**2)/(a**2-m**2))*P , Pi
    d = 1
    for a in range(1,m+1):
        d *= (2*a+1)/(2*a)  
    res = ms.log10(abs(Pi))+m/2.0*ms.log10(1-x**2)+1/2.0*ms.log10(d/(4*ms.pi))
    return res

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
    plt.title(r"Wave function with n = "+str(n)+r" m = "+str(m)+" and l = "+str(l))
    plt.show()

def Output_Harmonic(filename):
    '''
    测试球谐函数
    '''
    phi = ms.pi/5.0
    for theta in [ms.pi/1000,3*ms.pi/10,501*ms.pi/1000]:
        for L in [100,500,1000]:
            for M in [1,L/100.0,L/10.0,L-1]:
                try:
                    Y = New_GeneralizedLegrend(int(L),int(M),ms.cos(theta))
                except OverflowError as p:
                    print(L,M,theta)
                print("theta = "+str(theta)+" L = "+str(L)+" M = "+str(M)+\
                    ", Log_10(Y) = "+str(Y),file=filename)

def energy():
    '''
    能量函数，计算氢原子能量
    '''
    epsi = 10**(-6)# 数值精度
    d =60.0 # 数据范围
    n = 2 # 初始分隔
    # SP = d*y(d/2)
    S = d/n*(y(0)+y(d))
    while abs(S*8/3**9+1/18.0) > epsi:
        n *= 2 
        S = sum([y(d*i/(n-1)) for i in range(0,n)])/(n/d)
        print('n = '+str(n))
        print('Error = '+str(S*8/3**9+1/18.0))
    return S*8/3**9


def y(r):
    h = 0.001
    k = f(r)
    dy = (f(r+h)-f(r-h))/(2*h)
    ddy = (f(r+h)+f(r-h)-2*f(r))/h**2
    return k*(k*(1-r)-r*dy-r**2*ddy/2)

def f(r):
    return ms.exp(-r/3.0)*r*(6-r)
if __name__ == '__main__':
    # print(LaguerreL(3,20,10**(-3)))
    # print(Harmonic(100,99,ms.pi/1000,ms.pi/5))
    # print(New_GeneralizedLegrend(1000,999,ms.cos(ms.pi/1000)))
    # draw(2,1,1,5)
    # test(2,1)
    # f = open(r"./data/Harmonic.txt","w")
    # Output_Harmonic(f)
    # f.close()
    print(energy())
