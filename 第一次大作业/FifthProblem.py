'''
正定带状矩阵进行Chelosky分解，要求时空复杂度做到最小
'''
# 矩阵半带宽为2

n = 100 # 矩阵大小

# 利用三个向量存储矩阵
a = [[6 for _ in range(n-1)],[4 for _ in range(n)],[1 for _ in range(n-1)]]
a[0].append(5)
a[0].insert(0,5)

def matrix(i,j,mat): # 定义取矩阵算法
    m = i- j
    if m >= 0 and m < 3:
        return mat[m][i]
    else :
        raise Exception("Matrix error m ="+str(m))

