'''
正定带状矩阵进行Chelosky分解，要求时空复杂度做到最小
'''
# 矩阵半带宽为2
from numpy import *
from math import *

n = 100 # 矩阵大小
m = 2 # 带宽

# 利用三个向量存储矩阵

def solve(a,b):
    for i in range(n): # 下三角前代过程
        b[i] -= sum([a[i-j][j]*b[j] for j in range(max(0,i-m),i)])
        b[i] /= a[0][i]
    for x in range(n-1,-1,-1):
        b[i] -= sum([a[j-i][i]*b[j] for j in range(i+1,min(n,i+m))])
        b[i] /= a[0][i]
    return b

def validation(n, b):
  A = diag([5]+([6]*(n-2))+[5])
  A += diag([4]*(n-1), k=1)
  A += diag([4]*(n-1), k=-1)
  A += diag([1]*(n-2), k=2)
  A += diag([1]*(n-2), k=-2)
  return A.dot(asarray(b))

if __name__ =="__main__":
    a = [[5]+[6 for _ in range(n-2)]+[5],[4 for _ in range(n-1)],[1 for _ in range(n-2)]]
    b = [60]+([120]*(n-2))+[60]

    
    a[0][0]= sqrt(a[0][0])
    for i in range(1,n): # 行
        for j in range(max(0,i-m),i): # 列
            a[i-j][j] -= sum([a[i-k][k]*a[j-k][k] for k in range(max(0,i-m),j)])# 计算列中元素
            a[i-j][j] /= a[0][j]
        a[0][i] -= sum([a[i-k][k]**2 for k in range(max(0,i-m),i)])# 先计算对角元素
        try :
            a[0][i] = sqrt(a[0][i])
        except ValueError as p : # 防止产生错误
            print (p)
    print(a[0])
    print(a[1])
    print(a[2])
    b=solve(a,b)
    # print(b)
    print(len(validation(n,b)))