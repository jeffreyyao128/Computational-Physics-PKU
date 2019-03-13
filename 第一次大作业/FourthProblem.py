#coding=utf-8
'''
题目的要求是用三次样条插值内插出飞机的曲线
程序分为求解和画图两个部分
题目所给条件为自然边界条件
'''

x0 = [0,3,5,7,9,11,12,13,14,15] # 函数的定义域为[0,15]
y0=[0,1.2,1.7,2.0,2.1,2.0,1.8,1.2,1.0,1.6] 
n =len(x)

def firstorder(x,y):
    '''
    返回插值中n个点的一阶导数m_i
    采用三转角方程求解一阶导数
    '''

    pass

def f(x,x0,y0,m):
    '''
    返回三次样条内插函数值
    '''
    if x < x0[0] and x>x0[n-1]:
        raise ValueError('You x value is out of range!!')
    for i in range(1,n):
        if x < x0[i] and x > x0[i-1]:
            h = x[i]-x[i-1]
            x1 = x - x[i-1]
            x2 = x[i] -x
            return ((h+2*x1)*x2**2*y0[i-1]/h**3+\
            (h+2*x2)*x1**2*y0[i]/h**3+\
            x1*x2**2*m[i-1]/h**2-\
            x2*x1**2*m[i]/h**2)

if __name__=='__main__':
    m = firstorder(x0,y0)

