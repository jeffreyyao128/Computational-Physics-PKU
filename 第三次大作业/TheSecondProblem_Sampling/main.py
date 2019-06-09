# encoding=utf-8
'''
第三次大作业第二题
半经典激光电离隧穿问题
'''

import random
import numpy as np
import matplotlib.pyplot as plt

class Field:
    '''
    电场类
    '''
    def __init__(self,o=0.057,A0=1.325,type0='lin'):
        self.omega = o
        self.A0 = A0
        self.type0 = type0
    def E(self,t):
        o = self.omega
        if self.type0=='lin':
            return [-self.A0*o*(np.cos(o*t/8)**2*np.cos(o*t)-.25*np.cos(o*t/8)*np.sin(o*t/8)*np.sin(o*t)), 0]
        elif self.type0 == 'elip':
            return
    def E_amp(self,t):
        return np.sqrt(sum([each**2 for each in self.E(t)]))

class sample:
    '''
    取样类
    '''
    omega = 0
    A0 = 0
    def __init__(self,t,o=0.057,A0=1.325):
        self.omega = o
        self.A0 = A0
        self.t = t
        self.F = Field(o,A0)
        self.E_t=np.sqrt(sum([each**2 for each in self.F.E(t)]))

    def gaussian(self):
        '''
        采用第二类舍选法采样一次
        输出采样电子速度v_h的模值
        '''
        repeat = True
        while repeat:
            s1 = random.uniform(0,1)
            s2 = random.uniform(0,1)
            x = - np.log(s1)
            if (x - self.E_t/2)**2<-self.E_t*np.log(s2):
                return x
            # print('again')


class electron:
    '''
    电子类
    '''
    q = [0,0]
    def __init__(self,v,t0,type0,o=0.057,A0=1.325):
        self.t0 = t0
        self.p = v
        self.F = Field(o,A0,type0)
        self.q = [-0.5*self.F.E(t0)[0]/self.F.E_amp(t0)** 2, -0.5*self.F.E(t0)[1]/self.F.E_amp(t0)**2]
    
    def __p_evolve(self,q,t):
        return [-self.F.E(t)[0]-q[0]/(q[0]**2+q[1]**2+0.2**2)**1.5, -self.F.E(t)[1]-q[1]/(q[0]**2+q[1]**2+0.2**2)**1.5]
    
    def __q_evolve(self,p,t):
        return p[:]
    def evolve(self,dt = 0.1,T=2*np.pi/0.057):
        '''
        利用RK4进行牛顿力学模拟
        '''
        n=2
        t_s = np.linspace(self.t0,2*T,int((2*T-self.t0)/dt))
        # x,y = [],[]
        for t in t_s:
            # RK4计算哈密顿方程
            kq1 = self.__q_evolve(self.p,t)
            kp1 = self.__p_evolve(self.q,t)
            kq2 = self.__q_evolve([self.p[i]+.5*dt*kp1[i] for i in range(n)],t+0.5*dt)
            kp2 = self.__p_evolve([self.q[i]+.5*dt*kq1[i] for i in range(n)],t+0.5*dt)
            kq3 = self.__q_evolve([self.p[i]+.5*dt*kp2[i] for i in range(n)],t+0.5*dt)
            kp3 = self.__p_evolve([self.q[i]+.5*dt*kq2[i] for i in range(n)],t+0.5*dt)
            kq4 = self.__q_evolve([self.p[i]+dt*kp3[i] for i in range(n)],t+dt)
            kp4 = self.__p_evolve([self.q[i]+dt*kq3[i] for i in range(n)],t+dt)
            # 更新p,q
            self.p = [self.p[i]+dt/6*(kp1[i]+2*kp2[i]+2*kp3[i]+kp4[i]) for i in range(n)]
            self.q = [self.q[i]+dt/6*(kq1[i]+2*kq2[i]+2*kq3[i]+kq4[i]) for i in range(n)]
            # x.append(self.q[0])
            # y.append(self.q[1])
        #下面是测试代码
        # plt.plot(x,y)
        # plt.show()
        pf = np.sqrt(sum([each**2 for each in self.p]))
        rf = np.sqrt(sum([each**2 for each in self.q]))
        p_inf = np.sqrt(pf**2-2/rf)
        L = self.q[0]*self.p[1]-self.q[1]*self.p[0] # z右手方向为正
        a = [self.p[1]*L-self.q[0]/rf,-self.p[0]*L-self.q[1]/rf]
        return [p_inf*(p_inf*(-L*a[1])-a[0])/(1+p_inf**2*L**2), p_inf*(p_inf*(L*a[0])-a[1])/(1+p_inf**2*L**2)]



def test_gaussian(N):
    test = sample(1)
    test.E_t = 1
    y = [0 for i in range(300)]
    for _ in range(N):
        p = test.gaussian()
        if p > 3:
            continue
        else:
            y[abs(int(p*100))]+=1
    y = np.array(y)/y[0]
    x = np.linspace(0,3,300)
    y0 = np.exp(-x**2/test.E_t)
    plt.plot(x,y,label='sampling')
    plt.plot(x,y0,label='Calc')
    plt.show()

def lin_cal(dp=0.02,dt=1,T = 2*np.pi/0.057,N=10**5):
    t_s = np.linspace(-2*T,2*T,int(4*T/dt))
    type0 = 'lin'
    Elec = Field(type0=type0)
    px = [0 for i in range(int(3/dp))]
    py = [0 for i in range(int(3/dp))]
    data = [[0 for _ in range(int(3/dp))] for _ in range(int(3/dp))]
    for t in t_s:
        magnitude = Elec.E_amp(t)
        if magnitude<0.03:
            continue
        print("t = "+str(t))
        n = N*np.exp(-2*magnitude/3)/(100*magnitude)**1.5
        print("n = "+str(n))
        for i in range(int(n)): #进行n次采样
            S = sample(t)
            v_mag = S.gaussian()
            v = [v_mag*Elec.E(t)[1]/magnitude,-v_mag*Elec.E(t)[0]/magnitude]
            e = electron(v,t,type0)
            p_res = e.evolve(dt = 3*dt,T=T)
            if p_res[0]<1.5 and p_res[0]>-1.5:
                px[int((p_res[0]+1.5)/dp)]+=1
            if p_res[1] < 1.5 and p_res[1] > -1.5:
                py[int((p_res[1]+1.5)/dp)]+=1
            if p_res[0] < 1.5 and p_res[0] > -1.5 and p_res[1] < 1.5 and p_res[1] > -1.5:
                data[int((p_res[0]+1.5)/dp)][int((p_res[1]+1.5)/dp)]+=1
            if i%100==0:
                print(i)   
    x = np.linspace(-1.5,1.5,int(3/ dp))
    plt.figure(1)
    plt.plot(x,px)
    plt.plot(x,py)
    plt.figure(2)
    fig, ax =plt.subplots()
    im = ax.imshow(data)
    plt.colorbar(im)
    plt.show()
                     


if __name__ == '__main__':
    # test_gaussian(10**6)
    # e = electron([0,1],0,'lin')
    # e.evolve()
    lin_cal(N=10000,dt=0.5)

