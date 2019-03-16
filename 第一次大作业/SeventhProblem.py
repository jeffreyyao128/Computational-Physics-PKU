#coding=utf-8
'''
计物大作业第七题
求解定态氢原子波函数
'''

import math as ms

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

if __name__ == '__main__':
    # print(LaguerreL(30,40,100))
    # print(Harmonic(2,1,ms.pi/3,-ms.pi/3))
    # print(GeneralizedLegrend(2,1,0.5))
    