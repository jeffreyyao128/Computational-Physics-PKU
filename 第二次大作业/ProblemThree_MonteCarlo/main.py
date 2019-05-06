#encoding=utf-8
'''
第三题，蒙特卡洛问题
'''

import numpy as np
from time import time

def pseudo_rand(n=1,x0 = int(time())):
    '''
    线性同余伪随机数产生
    '''
    m = 2**31-1
    a = 1103515245
    c = 12345
    res = []
    for _ in range(n):
        x0 = (a*x0+c)%m
        res.append(x0)
    res = [each/float(m) for each in res]
    return res

def walk(x):
    if x>=0.75 and x<1:
        return (1,0) # 向右
    elif x<0.75 and x>=0.5:
        return (0,1) # 向上
    elif x<0.5 and x>=0.25:
        return (-1,0) # 向左
    elif x>=0 and x<0.25:
        return (0,-1) # 向下
    else :
        raise('Range Error')

def rand_walk(n=50):
    '''
    第二问
    '''
    a = pseudo_rand(n)
    s = [walk(x) for x in a]
    x = y =0
    res = []# 保存路径
    for w in s:
        x+=w[0]
        y+=w[1]
        res.append([x,y])
    return res

def Third_rand_walk(n=100):
    '''
    第三问
    '''
    # a = [pseudo_rand()[0] for _ in range(n)]
    # a = pseudo_rand(n)
    a = [np.random.random_sample() for _ in range(n)]
    sequ = [2*np.pi*each for each in a]
    x=y=0
    # res = []
    for each in sequ:
        x+=np.cos(each)
        y+=np.sin(each)
        # R = np.sqrt(x**2+y**2)
        # res.append([x,y,R])
    # return res
    return np.sqrt(x**2+y**2)

import matplotlib.pyplot as plt

def Third_plot(m=1000):
    '''
    第三问画图
    '''
    sqrtN =[]
    R =[]
    for x in range(11):
        sqrtN.append(np.sqrt(10*x))
        # sqrtN.append(10*x)
        avr = sum([Third_rand_walk(10*x) for _ in range(m)])/m # 取样m次做平均
        R.append(avr)
    # print(R)
    x = np.array(sqrtN)
    y = np.array(R)
    plt.scatter(x,y,alpha=0.6)
    plt.xlabel(r'$\sqrt{N}$')
    plt.ylabel(r'$\bar{R}$')
    plt.title(r'Relation between $\bar{R}$~$\sqrt{N}$')
    # 拟合
    a = sum([x[i]*y[i] for i in range(len(x))])/sum([each**2 for each in x])
    r = sum([x[i]*y[i] for i in range(len(x))]) / \
        np.sqrt(sum([each**2 for each in x])*sum([each**2 for each in y]))
    x1 = np.linspace(0,10)
    y1 = a*x1
    plt.plot(x1,y1)
    plt.show()
    print("a = "+str(a))
    print("r = "+str(r))





if __name__ == '__main__':
    # print(pseudo_rand(5))
    # f = open(r'./res.txt','w')
    # for each in rand_walk(50):
    #     x = each[0]
    #     y = each[1]
    #     print("x ="+str(x)+",y ="+str(y)+",R ="+str(np.sqrt(x**2+y**2)),file=f)
    # f.close()
    # f = open(r'./res_2.txt','w')
    # for each in Third_rand_walk():
    #     x = each[0]
    #     y = each[1]
    #     R = each[2]
    #     print("x ="+str(x)+",y ="+str(y)+",R ="+str(R),file=f)
    # f.close()
    Third_plot()

    

