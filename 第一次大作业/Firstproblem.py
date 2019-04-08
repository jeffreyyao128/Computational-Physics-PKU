'''
问题描述：编写两个程序，实现任意十进制数与其规范化的二进制双精度浮点数编码间的
相互转换，并列举足够算例说明你程序的正确性。
'''

import numpy as np # 用于判断符号位
import math as ms # 用于计算对数和向下取整

# 首先实现二进制转十进制
def bi2deci(a):
    '''
    实现double型二进制数转十进制 
    输入为a[64]数据类型为<class ='list'>
    输出为float(python默认为double)
    '''
    sign = -1 if a[0]==1 else 1 #符号位
    E = 0# 指数位
    for i in range(1,12):
        E += 2**(11-i)*a[i]
    # print(E)
    T = 1# 尾数位
    for i in range(12,64):
        T += a[i]*2**(11-i)
    # print(T)
    x = sign * T * 2**(E - 1023)
    if x < 2**(-1022):
        return 0
    return x


# 十进制转二进制

def deci2bi(a):
    '''
    实现double(float)型十进制数转二进制
    返回res[64]
    a[0]为符号,a[1:12]为指数位,a[12:64]为小数位
    '''
    S = np.sign(a)
    res = [0 for i in range(64)]
    if S ==0:
        return res
    S = 0 if np.sign(a)==1 else 1
    res[0] = int(S)
    a = abs(a)
    E = ms.floor(ms.log2(a))
    if a >=2:
        # 符号位> 1023
        while a>=2 :
            a/=2
            # E+=1
    elif a<1:
        while a <1:
            a*=2
            # E-=1
    # 此时 a \in (1,2]
    a -= 1 #去掉首位的1
    # print("E = "+str(E))
    E +=1023 #加上偏移量
    for k in range(11,0,-1): #计算指数位
        res[k] =int(E %2.0)
        E /=2.0
        if E==0:
            break
    p = 1
    while p<=52:
        a *=2
        h = 1 if a>=1.0 else 0
        res[11+p] = h
        a -=h
        p+= 1
    return res

if __name__=='__main__':
    p=(deci2bi(0.5))
    for i in p:
        print(i,end='')
    print()
    print(bi2deci(p))