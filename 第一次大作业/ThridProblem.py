#encoding=utf-8
'''
第三题：求解耦合的衰变微分方程
利用差分法递推求解
'''
import matplotlib.pyplot as plt
import numpy as np # 这两个包仅用作画图

def decay(ta=1,tb=1,dt=0.01,tmax=1):
    '''
    衰变函数
    返回<class='tuple'> (NA,NB)
    '''
    n = 0
    NA=[1]
    NB=[1]
    while n*dt < tmax :
        NA.append(NA[n]*(1-dt/ta))  
        NB.append(NB[n]*(1-dt/tb)+NA[n]*(dt/ta))
        n+=1
    return (NA,NB)

def draw(NA,NB,ta=1,tb=1,dt=0.01,tmax=1):
    '''
    作图函数
    '''
    # plt.rc('text',usetex=True)
    t = np.linspace(0,tmax,int(tmax/dt)+1)
    na = np.array(NA)
    nb = np.array(NB)
    test = np.exp(-t/ta)
    # plot N_A
    plt.subplot(2,1,1)
    plt.title(r'Decay with $\tau_A$ = '+str(ta)+r' and $\tau_B  = $'+ str(tb))
    plt.xlabel('time')
    plt.ylabel('$N_a$')
    plt.plot(t,na,'r-',label='Calculative')
    plt.plot(t,test,'--',label='Theoretical')
    plt.legend()
    # plot N_B
    plt.subplot(2,1,2)
    plt.plot(t,nb,'r-')
    plt.xlabel('time')
    plt.ylabel('$N_b$')
    plt.show()

def dif():
    '''
    计算相对误差函数并作图
    '''
    tmax = 10
    ta , tb = 1, 10

    plt.subplot(3,1,1)
    plt.title('Absolute errors for different $\Delta t$')
    dt = 0.2
    NA ,NB = decay(ta,tb,dt,tmax) 
    t = np.linspace(0,tmax,int(tmax/dt)+1)
    na = np.array(NA)
    nb = np.array(NB)
    test = np.exp(-t/ta)
    draw = np.abs((test-na))
    plt.plot(t,draw,'r--',label='$\Delta t = 0.2$')
    plt.legend()

    plt.subplot(3,1,2)
    dt = 0.1
    NA ,NB = decay(ta,tb,dt,tmax) 
    t = np.linspace(0,tmax,int(tmax/dt)+1)
    na = np.array(NA)
    nb = np.array(NB)
    test = np.exp(-t/ta)
    draw = np.abs((test-na))
    plt.plot(t,draw,'g--',label='$\Delta t = 0.1$')
    plt.legend()
    
    plt.subplot(3,1,3)
    dt = 0.05
    NA ,NB = decay(ta,tb,dt,tmax) 
    t = np.linspace(0,tmax,int(tmax/dt)+1)
    na = np.array(NA)
    nb = np.array(NB)
    test = np.exp(-t/ta)
    draw = np.abs((test-na))
    plt.plot(t,draw,'b--',label='$\Delta t = 0.05$')
    plt.legend()
    
    plt.show()

if __name__=='__main__':
    # ta = 1
    # tb = 10
    # dt = 0.2
    # tmax = 100
    #NA,NB = decay(ta=ta,tb=tb,dt=dt,tmax=tmax)
    # draw(NA,NB,ta=1,tb=0.1,dt=0.01,tmax=100)
    dif()



