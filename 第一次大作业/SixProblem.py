'''
利用共轭梯度法迭代求解稀疏矩阵的线性方程组
由于迭代法实际上只需要进行矩阵乘法，因此只用定义乘法函数
'''
from scipy.sparse import csr_matrix # 使用了scipy库中的行压缩稀疏矩阵
import numpy as np # 用作矩阵乘法
import gc #用于内存管理
import matplotlib.pyplot as plt # 用作画图


e =  10**(-6) #最终精度
n = 100 #矩阵大小

def CGM(e,b,A,x=[0 for _ in range(n)]):
    '''Input error bound e and ndarray b and sparse matrix A
    return solution ndarray x
    '''
    x = np.array(x)
    r = b - A*x # 计算初始残差
    p = r # 计算初始方向矢量
    epsilon = np.dot(r.T,r)
    k = 0
    while epsilon >= e :
        delt = A*p
        a = np.dot(r,p)/np.dot(r,delt) # a是标量
        x = x + a*p
        r = r - a*delt
        beta = - np.dot(r,delt)/np.dot(p,delt)
        p = r + beta*p
        k += 1
        epsilon = np.dot(r.T,r)
        print(k,"error : ",epsilon)
    return x

def draw(res):
    '''
    绘图函数
    '''
    x = np.linspace(0,n-1,n)
    y = np.array(res)
    plt.plot(x,y,'go',label='Solution')
    plt.xlabel('i')
    plt.ylabel('x[i]')
    plt.legend()
    plt.show()

if __name__ == "__main__" :
    b = np.array([2.5]+[1.5]*int(n/2-2)+[1.0]*2+[1.5]*int(n/2-2)+[2.5],dtype = float)
    a1 = np.diag([3]*n,k=0).astype(np.int8)
    a2 =a1 + np.diag([-1]*(n-1),k=1).astype(np.int8)
    del a1
    gc.collect()
    a3 =a2 + np.diag([-1]*(n-1),k=-1).astype(np.int8)
    del a2
    gc.collect()
    a4 = np.rot90(np.diag([0.5]*int(n/2-1)+[0,0]+[0.5]*int(n/2-1)),k=0).astype(np.float16)
    A = a3 +a4
    A = csr_matrix(A)

    x = CGM(e,b,A)
    draw(x)
    # print(x)
    # print(A*x)