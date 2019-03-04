'''
问题描述：编写两个程序，实现任意十进制数与其规范化的二进制双精度浮点数编码间的
相互转换，并列举足够算例说明你程序的正确性。
'''

# 首先实现二进制转十进制

def bintodeci(str):
    (pre_part,_,post_part)=str.partition(".") # 将数从小数点分成三个部分
    res = 0
    if _ == '.':        # 输入中存在小数点
        for index in range(len(pre_part)): # 整数部分
            a = int(pre_part[-index-1]) # 倒着排
            if a in [1,0]:
                res += a * (2 ** index)
            else :
                return "wrong input!"
        for index in range(len(post_part)): # 小数部分
            a = int(post_part[index])
            if a in [1,0]:
                res += a * (2 ** (-index -1))
            else :
                return "wrong input!"
        return res
    elif _ == '' : # 不存在小数点
        for index in range(len(pre_part)): # 整数部分
            a = int(pre_part[-index-1]) # 倒着排
            if a in [1,0]:
                res += a * (2 ** index)
            else :
                return "wrong input!"
        return res


# 十进制转二进制

def decitobin(string):
    p=False
    while p == False:
        digits = input ("How many digits of fraction part do you want to display? :\n")    # 确定显示位数
        try:
            digits = int(digits)
        except ValueError as e:
            print ("Warrning"+e)
            continue
        if digits > 0 :
            p = True
    (pre_part,_,post_part)=string.partition(".")
    integer=""
    decimal=""
    if _ == '.': #输入中存在小数点
        # 整数部分不停除二取余
        pre=int(pre_part)
        while pre != 0:
            re = str(pre%2)
            integer =integer+ re
            pre = pre//2
        #小数部分不停乘二取整
        post=float("0"+_+post_part)
        while post != 0 and len(decimal) < digits :
            re = str(int(post*2))
            decimal = decimal+ re
            post = (post*2)%1 # 取小数位
        if len(decimal)<digits: # 补齐小数位
            decimal = decimal + '0'*(digits-len(decimal))
    elif _ == "": # 输入中不存在小数点
        # 整数部分不停除二取余
        pre=int(pre_part)
        while pre != 0:
            re = str(pre%2)
            integer =integer+ re
            pre = pre//2
    result = integer+_+decimal
    return result

if __name__=='__main__':
    end =False
    while end != True:
        a=input("Binary or Decimal?: Bin/Deci \n")
        if a.lower() == 'bin':
            num =input ("please input your number:")
            print(bintodeci(num))
            print("=================")
        elif a.lower() == 'deci':
            num =input("please input your number:")
            print(decitobin(num))
            print("=================")
        elif a.lower() == "exit":
            break
        else :
            print("please try again!")
            continue