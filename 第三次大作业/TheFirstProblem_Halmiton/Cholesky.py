'''
Cholesky 方法
'''
import numpy as np

def cholesky(a):
    """
    Cholesky求解方法,返回矩阵L
    return a
    """
    m=len(a)-1
    n = len(a[0])
    # print(a)
    for k in range(n):
        #先解对角元素
        a[0][k] -= sum(
            [a[k - j][j]**2 * a[0][j] for j in range(max(0, k - m), k)])
        for l in range(1, min(n - k, m + 1)):
            a[l][k] -= sum([
                a[k - j][j] * a[0][j] * a[l + k - j][j]
                for j in range(max(0, l + k - m), k)
            ])
            a[l][k] /= a[0][k]
    return a

def solve(a,b):
    '''
    解方程方法
    返回解向量
    '''
    m = len(a) - 1
    n = len(a[0])
    # print(b)
    for k in range(n):
        b[k] -= sum([a[k-j][j]*b[j] for j in range(max(0, k-m), k)])
    for k in range(n):
        b[k] /= a[0][k]
    for k in range(n):
        ki = n-1-k
        b[ki] -= sum([a[j][ki]*b[ki+j] for j in range(1, min(n-ki, m+1))])
    return b

if __name__=="__main__":
    '''
    测试部分
    '''
    x0=.48
    xmax =2000
    N =200
    Dx = 2*xmax/(2*N-1)
    A = [[1/Dx**2-1/np.sqrt((i*Dx)**2+2)+x0 for i in np.linspace(-xmax,xmax, N)], [-1/(2*Dx**2) for i in range(N-1)]]  # 已经加上位移
    D1 = np.diag([1/Dx**2-1/np.sqrt((i*Dx)**2+2)+x0
                  for i in np.linspace(-xmax, xmax, N)])
    D1 += np.diag([-1/(2*Dx**2) for i in range(N-1)], k=1) + \
        np.diag([-1/(2*Dx**2) for i in range(N-1)], k=-1)
    L = cholesky(A)
    # print(L)
    b = np.array([.5 for i in range(N)])
    v = np.array(solve(L,b.copy()))
    a = D1@v
    print(a)
    print(b)
    print(np.linalg.norm(a-b))
