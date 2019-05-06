# encoding=utf-8
'''
第一题：解决Paul Trap 中的mathium equation的本征值问题
注意：题目最后一问要求求解本征矢量
'''

import numpy as np
import scipy.linalg as lin# 仅用于检验算法的正确性
import matplotlib.pyplot as plt# 用于作图
from scipy import interpolate

def construct(M,q):
    '''
    构造H矩阵,输入参数M = 
    '''
    N = float(2*M+1)
    V = np.diag([2*q*np.cos(2*2*np.pi*k/(N)) for k in range(2*M+1)]) # Diagonal term
    T = np.array([\
        [((-1)**(j-k)*np.cos((j-k)*np.pi/N)/(2*np.sin(np.pi*(j-k)/N)**2) if j!=k else M*(M+1)/3.0) for j in range(2*M+1)]\
        for k in range(2*M+1)])
    return (V+T)

def norm(x):
    return np.sqrt(np.sum(np.square(x)))

def veci(func):
    return np.vectorize(func,otypes=[int])

@veci
def filtrate_i(x,e=10**(-15)):
    '''
    过滤函数，用于过滤掉机器精度以下的数据
    '''
    return x if abs(x)>e else 0


def vecf(func):
    return np.vectorize(func, otypes=[float])


@vecf
def filtrate_f(x, e=10**(-15)):
    '''
    过滤函数，用于过滤掉机器精度以下的数据
    '''
    return x if abs(x) > e else 0


def Householder(A):
    '''
    Householder变换，返回Hessengerg矩阵和变换矩阵(A,G)
    输入A:矩阵
    '''
    n = A.shape[0]
    P = np.identity(n)
    for i in range(n-2):
        # 构建反射矩阵
        x = A[i+1:,i] # 取矩阵列共(n-i-1)维
        e = np.array([1]+[0 for x in range(n-i-2)]) # 构造单位矢量
        v = np.sign(x[0])*norm(x)*e+x # 构造反射矢量
        v = v/norm(v) #归一化 
        G = np.array([[2*x*y for x in v]for y in v])
        H = lin.block_diag(np.identity(i+1),np.identity(n-i-1)-G)# 构建反射矩阵
        A = H@A@H
        P = H@P
    # A = filtrate(A)
    return A,P

def decompose(A):
    '''
    利用Givens变换分解，A为上海森堡矩阵
    返回(Q,R)
    '''
    n = A.shape[0]
    G = np.identity(n) # Givens 矩阵
    for i in range(n-1): # 共(n-1)个非对角元需要消去
        cos = A[i,i]/np.sqrt(A[i,i]**2+A[i+1,i]**2)
        sin = A[i+1,i]/np.sqrt(A[i,i]**2+A[i+1,i]**2)
        G1 = lin.block_diag(np.identity(i),np.array([[cos,sin],[-sin,cos]]),np.identity(n-i-2))
        A = G1@A
        G = G@G1.T
    return G , A

def is_schur(A,e=10**(-14)):
    '''
    函数判断是否为实舒尔型
    默认A为对称矩阵,精度为e
    '''
    # 这个估摸着用不上
    n = A.shape[0]
    for i in range(n-1):
        if A[i+1,i]>e and A[i+2,i+1]>e: return False
    return True


def QR(A):
    '''
    QR算法，输入Hessengberg矩阵A
    返回特征值和特征向量
    '''
    n = A.shape[0]
    i = 0
    N = 10**5
    e = 10**(-10)
    # G = np.identity(n)
    k = n-1
    res = [] # 储存本征矢量
    Q0 = np.identity(n)
    while k > 0:
        if abs(A[k,k-1])<e: # 本征值必然是实数
            res.append(A[k,k]) # 本征值代入
            A = A[:k,:k] # 取子阵
            k = k-1
            # print(k)
        s = A[k,k] # 位移因子
        A=A-s*np.identity(k+1)
        Q, R = decompose(A)
        A = R@Q + s*np.identity(k+1)
        Q0 = Q0@lin.block_diag(Q,np.identity(n-Q.shape[0]))
        i+=1
        if i > N: 
            raise("Too many Loops!")
    res.append(A[0,0])
    re=[]
    for i in range(n):
        re.append((res[n-1-i],Q0[:,i]))
    re = sorted(re, key=lambda s: s[0])
    return re

def solve(M,q):
    A = construct(M,q)
    G = Householder(A)[0]

    return [x[0] for x in QR(G)]

def draw(M):
    '''
    第二问绘图函数
    '''
    x = np.linspace(0,20,21)
    res = []
    for q in x:
        print(q)
        res.append(solve(M,q)[:11])
    y = np.array(res).T
    for i in range(11):
        plt.plot(x,y[i])
    plt.xlabel('q')
    plt.ylabel('Eigenvalues')
    plt.show()

def Even_vect(v):
    '''
    判断本征矢奇偶性
    '''
    n = v.shape[0]
    e = 10**(-2)
    if abs(v[1]-v[n-1])/abs(v[1])<e:
        return True
    elif abs(v[1]+v[n-1])/abs(v[1])<e:
        return False
    else:
        print(v)
        raise("Type error!")


def draw_third(M):
    '''
    第三题画图函数
    '''
    x = np.linspace(0,20,21)
    res = []
    for q in x:
        A = construct(M,q)
        G , H = Householder(A)
        s = QR(G)
        c =[]
        if q == 0:
            i = 0 
            while i < len(s):
                if abs(s[i][0]-s[i+1][0])<10**(-10): # 发生简并
                    c.append(s[i][0])
                    i+=1
                if len(c)==5:
                    break
                i+=1
        else:
            for a in s:
                if Even_vect(H.T@a[1]):
                    c.append(a[0])
                if len(c)==5:
                    break
        res.append(c[:5]) # 加入前五个
        print(q)
    y = np.array(res).T
    print(x)
    print(y)
    for i in range(5):
        plt.plot(x,y[i])
    plt.xlabel('q')
    plt.ylabel('Eigenvalues')
    plt.show()
        
def draw_Fourth():
    A = construct(50, 10)
    G, H = Householder(A)
    s = QR(G)
    c = []
    p =[]
    for a in s:
        if Even_vect(H.T@a[1]):
            c.append(a[0])
            p.append(a[1])
        if len(c) == 5:
            break
    # print(c)
    x = np.linspace(0,2*np.pi,2*50+1)[:10]
    for each in p:
        y = each[:10]
        tck, u = interpolate.splprep([x, y], s=0)
        xnew, ynew = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
        plt.plot(xnew,ynew)
    plt.show()

if __name__=="__main__":
    # print(sorted(QR(G)))    
    # print(sorted(lin.eig(G)[0]))
    # Q,R = decompose(G)
    # print(lin.eig(R@Q)[0])
    # Q,R = decompose(G)
    # print(Q@Q.T)
    #第二问
    # draw(50)
    #第三问
    # draw_third(5)
    # A = construct(5,5)
    # G ,H = Householder(A)
    # s = QR(G)
    # for x in s:
    #     print(lin.norm(G@x[1]/x[0]-x[1]))
    #     print(lin.norm(A@(H.T@x[1])/x[0]-H.T@x[1]))
    #     print(H.T@x[1])
    # 第四问
    draw_Fourth()
    
