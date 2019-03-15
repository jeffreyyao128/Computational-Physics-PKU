#coding=utf-8
'''
题目的要求是用三次样条插值内插出飞机的曲线
程序分为求解和画图两个部分
题目所给条件为自然边界条件
'''

import matplotlib.pyplot as plt
import numpy as np

x0 = [0,3,5,7,9,11,12,13,14,15] # 函数的定义域为[0,15]
y0=[0,1.2,1.7,2.0,2.1,2.0,1.8,1.2,1.0,1.6] 
n =len(x0)

def d(i):
    return (y0[i]-y0[i+1])/(x0[i]-x0[i+1])

def h(i):
    return (x0[i+1]-x0[i])


def firstorder(x,y):
    '''
    返回插值中n个点的一阶导数m_i
    采用三转角方程求解一阶导数
    '''
    # 建立矩阵
    mu=[h(i)/(h(i)+h(i+1)) for i in range(n-1)]
    lamb=[h(i+1)/(h(i)+h(i+1)) for i in range(n-1)]
    # 建立待解向量
    b = [3*d[0]]+[3*(lamb[i]*d(i)+mu(i)*d(i+1)) for i in range(n-1)]+[3*d[n-1]]

    
    return m

def f(x,m):
    '''
    返回三次样条内插函数值
    '''
    if x < x0[0] and x>x0[n-1]:
        raise ValueError('Your x value is out of range!!')
    for i in range(1,n):
        if x < x0[i] and x > x0[i-1]:
            h = x[i]-x[i-1]
            x1 = x - x[i-1]
            x2 = x[i] -x
            return ((h+2*x1)*x2**2*y0[i-1]/h**3+\
            (h+2*x2)*x1**2*y0[i]/h**3+\
            x1*x2**2*m[i-1]/h**2-\
            x2*x1**2*m[i]/h**2)

if __name__=='__main__':
    # m = firstorder(x0,y0)
    # 下面是画图部分
    x1 = np.linspace(0,15,100) # 在定义域中共100个点
    y1 = np.array([x **2 for x in np.nditer(x1)])
    
