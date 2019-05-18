'''
正定带状矩阵进行Chelosky分解，要求时空复杂度做到最小
'''
# 矩阵半带宽为2
import numpy as np  # 调用numpy仅在validation函数中用作检验
import math # 用于开根号计算
import matplotlib.pyplot as plt # 用作绘图

n = 1000 # 矩阵大小
m = 2 # 带宽

# 利用三个向量存储矩阵

def validation(n, b):
    '''
    验证解，调用numpy矩阵运算
    '''
    A = np.diag([5]+([6]*(n-2))+[5])
    A += np.diag([4]*(n-1), k=1)
    A += np.diag([4]*(n-1), k=-1)
    # A += np.diag([1]*(n-2), k=2)
    # A += np.diag([1]*(n-2), k=-2)
    return A.dot(np.asarray(b))

# def NeoCholesky(a,b):
#     '''
#     Cholesky 分解，返回解向量
#     '''
#     for i in range(n):
#         for j in range(max(0,i-m),i+1):
#             a[i-j][j] -= sum([a[i-k][k]*a[j-k][k]/a[0][k] for k in range(max(0,i-m),j)])
#     for i in range(n):
#         b[i] -= sum([a[i-j][j]*b[i]/a[0][j] for j in range(max(0,i-m),i)])
#     # 上三角回代过程
#     for i in range(n-1,0,-1):
#         try:
#             b[i] -= sum([a[j-i][i]*b[i] for j in range(i+1,min(n,i+m+1))])
#         except IndexError as p:
#             print(p)
#             print(" i = "+ str(i) )
#         b[i] = b[i] /a[0][i]
#     return a,b

def cholesky(a, b):
    """
    Cholesky求解方法,返回解向量和矩阵
    return a,b
    """
    # print(a)
    n= len(a[0])
    m =len(a)-1
    for k in range(n):
        #先解对角元素
        a[0][k] -= sum([a[k-j][j]**2*a[0][j] for j in range(max(0,k-m),k)])
        for l in range(1, min(n-k,m+1)):
            a[l][k] -= sum([a[k-j][j]*a[0][j]*a[l+k-j][j] for j in range(max(0,l+k-m),k)])
            a[l][k] /= a[0][k]
#回代解方程
    for k in range(n):
        b[k] -= sum([a[k-j][j]*b[j] for j in range(max(0,k-m),k)])
    for k in range(n):
        b[k] /= a[0][k]
    for k in range(n):
        ki = n-1-k
        b[ki] -= sum([a[j][ki]*b[ki+j] for j in range(1,min(n-ki,m+1))])
    return a,b

def draw(b,n):
    '''
    作图函数
    '''
    x = np.linspace(0,n-1,n)
    y = np.array(b)
    plt.plot(x,y,'r-')
    plt.xlabel('i')
    plt.ylabel("x[i]")
    plt.show()


if __name__ =="__main__":
    # a = [[5]+[6 for _ in range(n-2)]+[5],[4 for _ in range(n-1)]\
    #     # ,[1 for _ in range(n-2)]\
    #         ]
    # b = [60]+([120]*(n-2))+[63]
    # a1,b1=cholesky(a,b)
    # print(validation(n,b1))
    # draw(b1,n)
    A = [[3., 5., 3.], [1., 1.]]
    B = np.diag([3, 5, 3]) + np.diag([1, 1], k=1) + np.diag([1, 1], k=-1)
    b = np.array([1., 2., 3.])
    L,v = cholesky(A,b)
    # print(L)
    # print(np.linalg.cholesky(B))
    # v = np.array(solve(L, b))
    print(B@v)