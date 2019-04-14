# encoding=utf-8
'''
第一题：解决Paul Trap 中的mathium equation的本征值问题
注意：题目最后一问要求求解本征矢量
'''

import numpy as np
import scipy.linalg as lin

def construct(M,q):
    '''
    构造H矩阵
    '''
    V = np.diag([2*q*np.cos(2*np.pi*k/(2*M+1)) for k in range(2*M+1)]) # Diagonal term
    T = np.array([\
        [((-1)**(j-k)*np.cos((j-k)*np.pi/(2*M+1))/(2*(np.sin(np.pi*(j-k)/(2*M+1)))**2) if j!=k else M*(M+1)/3) for j in range(2*M+1)]\
        for k in range(2*M+1)])
    return (V+T)

def norm(x):
    return np.sqrt(np.sum(np.square(x)))

def vec(func):
    return np.vectorize(func,otypes=[int])

@vec
def filtrate(x,e=10**(-15)):
    '''
    过滤函数，用于过滤掉机器精度以下的数据
    '''
    return x if abs(x)>e else 0


def Householder(A):
    '''
    Householder变换，返回Hessengerg矩阵
    输入A:矩阵，n:维度
    '''
    n = A.shape[0]
    for i in range(n-2):
        # 构建反射矩阵
        x = A[i+1:,i] # 取矩阵列共(n-i-1)维
        e = np.array([1]+[0 for x in range(n-i-2)]) # 构造单位矢量
        v = np.sign(x[0])*norm(x)*e+x # 构造反射矢量
        v = v/norm(v) #归一化 
        G = np.array([[2*x*y for x in v]for y in v])
        H = lin.block_diag(np.identity(i+1),np.identity(n-i-1)-G)
        print(H)
        A = H@A@H
    A = filtrate(A)
    return A

def decompose(A):
    '''
    利用Givens变换分解，A为上海森堡矩阵
    返回(Q,R)
    '''
    pass

def QR(A):
    '''
    QR算法，返回特征值
    '''
    pass

if __name__=="__main__":
    A = construct(8,10)
    print(Householder(A))