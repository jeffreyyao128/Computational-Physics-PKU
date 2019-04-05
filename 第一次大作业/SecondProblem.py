#encoding=utf-8
'''
求解高斯求积的节点和求积系数
'''

import numpy as np #为了简单起见，这里调用了numpy的多项式类

def sqrtinte(p,ab=(0,1)):
    '''
    返回以\sqrt[x]为权函数的多项式p(x)的积分，上下限默认为(0,1)
    '''
    a = p.coef
    return sum([a[len(a)-i-1]/(i+3/2.) for i in range(0,len(a))])

def sqrtbuild(N):
    '''
    构造以\sqrt[x]为权函数的最高位N的正交多项式
    '''
    res = [0 for _ in range(N+2)] #把 n =-1存在了第一个
    res[0] = np.poly1d([0])
    res[1] = np.poly1d([1])
    for x in range(1,N+1):
        p1 = np.poly1d([1,0])*res[x]**2
        p2 = res[x]**2
        p3 = res[x-1]**2
        a = sqrtinte(p1)/sqrtinte(p2)
        b = sqrtinte(p2)/sqrtinte(p3) if x >1 else 0
        res[x+1] = np.poly1d([1,-a])*res[x] - b*res[x-1]
    return res[1:len(res)]

if __name__ == "__main__":
    # print(sqrtinte(np.poly1d([1,0])))
    # sqe = sqrtbuild(3)
    # print(sqrtinte(sqe[1]*sqe[0]))
    
        