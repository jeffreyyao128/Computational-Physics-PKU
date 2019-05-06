#encoding=utf-8
'''
第二次大作业第二题
利用鞍点法求解激光下的跃迁问题
'''
import math
import  numpy as np # 中间调用了poly1d多项式方法
import matplotlib.pyplot as plt # 用作绘图
from matplotlib.ticker import FuncFormatter # 仅用作绘图美化
import time

# 为了方便，将以下变量设为全局变量
omega = 45.5633525316/3200
A0 = np.sqrt(5*10**(13 - 16)/3.5094448314)/omega
Ip = 13.6/27.2
B = A0/2

def rand():
    return 4*np.random.random_sample()-2

def Newton_deflation(a,N=10**3,e=10**(-14)):
    '''
    牛顿法解方程
    输入a为多项式类<class 'numpy.poly1d'>
    N为最大迭代次数,e为迭代精度
    返回一个根和收缩多项式
    '''
    n = 0
    x = complex(rand(),rand())# 初始取值
    while n<N:
        d = np.poly1d((1,-x))
        q , b = a/d
        if abs(b[0])<e:
            # print(b)
            return x, q
        x-=b[0]/q(x)
    else:
        print(str(n)+" loops")
        return x, q
    
def find_roots(pz):
    '''
    找到方程的所有根
    a为多项式
    输入pz
    返回<class 'list'>
    '''
    D = 4/B*(np.sqrt(2*13.6/27.2) - pz*1j)
    a = np.poly1d((1, -2, 1, D, -1, 2, -1))  # 构造多项式
    res =[]
    for _ in range(6):
        x, a = Newton_deflation(a)
        res.append(x)
    t = [2*np.log(x)/(omega*1j) for x in res]
    for i in range(len(t)):
        t[i] = t[i] if t[i].imag > 0 else t[i].conjugate()
        while t[i]<0:
            t[i]+=4*np.pi/omega
        while t[i] > 4*np.pi/omega:
            t[i] -= 4*np.pi/omega
    return t,res

def check():
    '''
    检验函数，用于检验寻根的正确性
    '''
    pz=1
    D = 4/B*(np.sqrt(2*13.6/27.2) - pz*1j)
    a = np.poly1d((1, -2, 1, D, -1, 2, -1))  # 构造多项式
    t, x = find_roots(1)
    print(t)
    print([i for i in x])

def q(t,pz):
    return (B*np.sin(t)*(1-np.cos(t/2))+pz)

def dA(t):
    return 0.5*A0*np.sin(0.25*t)*(np.cos(t)*np.sin(0.25*t)+0.5*np.sin(t)*np.cos(0.25*t))


def S(tau,pz):
    res = - 2/3*B*pz + (0.75*B**2+pz**2)*tau
    res += 2*B*pz*(np.cos(tau/2)-np.cos(tau)+np.cos(1.5*tau)/3)
    res += B**2 * (-2*np.sin(tau/2)+np.sin(tau)/8+np.sin(1.5*tau)/3-0.375*np.sin(2*tau)+0.2*np.sin(2.5*tau)-np.sin(3*tau)/24)
    return (res/2 + Ip*tau)/omega 

def f(t,pz):
    '''
    返回被积函数
    '''
    return -(q(t,pz)*dA(t))/(q(t,pz)**2+2*Ip)**3*np.exp(1j*S(t,pz))

def SPM(pz):
    '''
    鞍点法求积积分
    '''
    t = find_roots(pz)[0]
    tau = [x*omega for x in t]
    M = 0
    for i in tau:
        M += np.exp(1j*S(i,pz))/(omega*q(i,pz)*dA(i))
    return (2*Ip)**(5/4.0)/np.sqrt(2)*M

def integ(pz,e=10**(-10),N=10**(6),a=0,b=4*np.pi):
    '''
    直接积分求解，输入pz
    输出积分值<class 'float'>
    '''
    n=2**8
    h = (b-a)/n
    T = 0.5*h*(f(a,pz)+2*sum([f(a+(i+0.5)*h,pz) for i in range(n)])+f(b,pz))
    while n<N:
        Tk=T
        T = T/2 + 0.5*h*sum([f(a+(i+0.5)*h, pz) for i in range(n)])
        n *= 2
        h /= 2
        error = abs(T-Tk)/abs(T)
        if error<e:
            return (4*T-Tk)/3
    else:
        print(str(n)+' loops!')
        return (4*T-Tk)/3

def draw_first():
    '''
    第一问画图函数
    '''
    p = np.linspace(0.01,2,100)
    E = 0.5*p**2
    M =[]
    for pz in p:
        M.append(abs(SPM(pz))**2)
        print(pz)
    M = np.array(M)
    f,a=plt.subplots(1,1)
    a.plot(E,M,'b')
    def formatnum(x, pos):
        return '$%.1f$x$10^{-5}$' % (x*10**(5))
    formatter = FuncFormatter(formatnum)
    a.yaxis.set_major_formatter(formatter)
    plt.xlabel("$E_k$")
    plt.ylabel("$|M_k|^2$")
    plt.title("SPM method")
    plt.show()

def draw_second():
    p = np.linspace(0.01, 2, 200)
    E = 0.5*p**2
    M = []
    for pz in p:
        # M.append(abs(integ(pz))**2)
        print(pz)
    M = np.array(M)
    f, a = plt.subplots(1, 1)
    a.plot(E, M, 'b')
    # def formatnum(x, pos):
    #     return '$%.1f$x$10^{-5}$' % (x*10**(5))
    # formatter = FuncFormatter(formatnum)
    # a.yaxis.set_major_formatter(formatter)
    plt.xlabel("$E_k$")
    plt.ylabel("$|M_k|^2$")
    plt.title("SPM method")
    plt.show()

if __name__=="__main__":
    # a = np.poly1d((1,0,1))
    # x,q = Newton_deflation(a)
    # print(x)
    # check()
    # print(abs(SPM(1))**2)
    # 第一问
    # for x in find_roots(1)[0]:
    #     print("Re = ",str(x.real),"Im = ",str(x.imag))
    # draw_first()
    # 第二问
    # print("M = ",abs(2**(3.5)*(2*Ip)**1.25/np.pi*integ(0.5))**2)
    # print(S(2,1))
    # draw_second()

    
