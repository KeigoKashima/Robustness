#!/usr/bin/env python3
from math import *
import numpy as np
from matplotlib import pyplot as plt

##ユニットの各パラメータ
M  = 1 #[kg] ユニットの質量
c  = 0.1    #粘性定数
k  = 0.1    #バネ定数
b  = 100.0/2*0.001  #[m] ユニットの幅/2
h  = 50.0/2 *0.001   #[m] ユニットの高さ/2
hd = 20.0  *0.001        #[m] バネダンパの作用点までの距離
I  = M*(b*b + h*h)/3  #１ユニットの重心まわりの慣性モーメント
g  = 9.8    ##[m/s2]重力加速度
x_max = 10.0*0.001    ##最大変位量
AllUnits = 20         #ユニット数

##[N]衝撃力
F = 10
Attacked_unit = 10

##時間
T  = 1.0    #[s] 全体の時間
TF = 0.02     #[s] 衝撃力がかかっている時間
dt = 0.01    #[s] 微小時間
step = int(T/dt) #ステップ数


##各変数の定義
##np.array(下から数えたユニットの番号，ステップ数)
x = np.zeros( ((AllUnits,step)))    ##ユニットiのx'座標: x[i][j]，
y = np.zeros ((AllUnits,step))    ##y'座標y[i][j]
xf = np.zeros( ((AllUnits,step))) ##x座標
yf = np.zeros ((AllUnits,step))    ##y座標

x_G  = np.zeros(step)
y_G  = np.zeros(step)
x_G0 = np.zeros(step) ##G0の重心x座標
x_G1 = np.zeros(step) ##G1
x_G2 = np.zeros(step) ##G2
y_G0 = np.zeros(step) ##G0の重心y座標
y_G1 = np.zeros(step) ##G1
y_G2 = np.zeros(step) ##G2
dx_G0 = np.zeros(step) ##G0の重心速度
dx_G1 = np.zeros(step) ##G1
dx_G2 = np.zeros(step) ##G2
ddx_G0 = np.zeros(step) ##G0の重心加速度
ddx_G1 = np.zeros(step) ##G1
ddx_G2 = np.zeros(step) ##G2

theta = np.zeros (step)##角度
w = np.zeros (step)    ##角速度
dw = np.zeros (step)   ##角加速度

M_G0 = np.zeros(step) ##G0の重心
M_G1 = np.zeros(step) ##G1
M_G2 = np.zeros(step) ##G2

I_G  = np.zeros(step)
I_G0 = np.zeros(step) ##G0の慣性モーメント
I_G1 = np.zeros(step) ##G1
I_G2 = np.zeros(step) ##G2

h_G0 = np.zeros(step) ##G0
h_G1 = np.zeros(step) ##重心
h_G2 = np.zeros(step) ##重心

##各点に働く力
f1 = np.zeros (step)   ##上の結合部に働く力
f2 = np.zeros (step)   ##下の結合部に働く力
N_G1 = np.zeros (step)
N_G2 = np.zeros (step)
n_G1 = np.zeros (step)
n_G2 = np.zeros (step)
FN   = np.zeros(step)
FH   = np.zeros(step)

##初期値
x[:,0] = -b
x_G0[0] = -b
x_G1[0] = -b
x_G2[0] = -b
y[:,0] += [i*2*h+h for i in range(AllUnits)]



def IF1F2(low_unit,high_unit,s):
    """
    引数
        low_unit：フリーユニット（下）
        high_unit：フリーユニット（上）
        s：ステップ数

    """
    num = high_unit - low_unit + 1

    M_G0[s] = M*num
    M_G1[s] = M*(AllUnits-high_unit-1)
    M_G2[s] = M*low_unit
    I_G0[s] = I_G0[s-1]
    I_G1[s] = I_G1[s-1]
    I_G2[s] = I_G2[s-1]
    ##左下まわりの全体の慣性モーメント
    I_G[s] = I_G0[s]+M_G0[s]*(x_G0[s]**2+y_G0[s]**2)\
             +I_G1[s]+M_G1[s]*(x_G1[s]**2+y_G1[s]**2)\
             +I_G2[s]+M_G2[s]*(x_G2[s]**2+y_G2[s]**2)
    y_G0[s] = y_G0[s-1]
    y_G1[s] = y_G1[s-1]
    y_G2[s] = y_G2[s-1]

    print("##dx_G1[s]",dx_G1[s],"x_G1[s]",x_G1[s],"y_G1[s]",y_G1[s],"M_G1[s]",M_G1[s],"I_G1[s]",I_G1[s])
    print("##dx_G0[s]",dx_G0[s],"x_G0[s]",x_G0[s],"y_G0[s]",y_G0[s],"M_G0[s]",M_G0[s],"I_G0[s]",I_G0[s])
    print("##dx_G2[s]",dx_G2[s],"x_G2[s]",x_G2[s],"y_G2[s]",y_G2[s],"M_G2[s]",M_G2[s],"I_G2[s]",I_G2[s])



    ##全体の重心位置
    x_G[s] = (M_G0[s]*x_G0[s]+M_G1[s]*x_G1[s]+M_G2[s]*x_G2[s])/(M*AllUnits)
    y_G[s] = h*AllUnits

    f1[s]= c*(dx_G0[s]-dx_G1[s]) + k*(x[high_unit,s]-x_G1[s])  #上
    f2[s]= c*(dx_G0[s]-dx_G2[s]) + k*(x[low_unit,s]-x_G2[s])   #下


def Inertia(low_unit,high_unit,s):
    """
    G0,G1,G2の重心座標，質量，慣性モーメント，垂直抗力，を求める関数
    引数
        low_unit：フリーユニット（下）
        high_unit：フリーユニット（上）
        s：ステップ数
        G0，G1，G2の重心，質量，慣性モーメントを求める
    """
    #間のユニットの数
    num = high_unit - low_unit + 1

    ##G0の重心の座標，質量
    x_G0[s] = sum(x[low_unit:high_unit+1, s])/num #x座標
    y_G0[s] = (y[high_unit+1,s]+y[low_unit+1,s])/2#y座標
    M_G0[s] = M*num
    h_G0[s] = h*num##2h*num/2
    I_G0[s] = 0
    for i in range(low_unit, high_unit+1):
        I_G0[s] += I + M*((x[i,s]-x_G0[s])**2+(y[i,s]-y_G0[s])**2)


    ##G1の重心の座標，質量
    x_G1[s] = sum(x[high_unit+1:, s])/(AllUnits-high_unit-1) #x座標
    y_G1[s] = sum(y[high_unit+1:, s])/(AllUnits-high_unit-1)#y座標
    M_G1[s] = M*(AllUnits-high_unit-1)
    h_G1[s] = h*(AllUnits-high_unit-1)
    I_G1[s] = M_G1[s]*(h_G1[s]**2+b**2)/3
    ##G2の重心の座標，質量
    x_G2[s] = sum(x[:low_unit, s])/low_unit #x座標
    y_G2[s] = sum(y[:low_unit, s])/low_unit#y座標
    M_G2[s] = M*low_unit
    h_G2[s] = h*low_unit
    I_G2[s] = M_G2[s]*(h_G2[s]**2+b**2)/3

    ##全体の重心位置
    x_G[s] = (M_G0[s]*x_G0[s]+M_G1[s]*x_G1[s]+M_G2[s]*x_G2[s])/(M*AllUnits)
    y_G[s] = h*AllUnits
    ##左下まわりの全体の慣性モーメント
    I_G[s] = I_G0[s]+M_G0[s]*(x_G0[s]**2+y_G0[s]**2)\
             +I_G1[s]+M_G1[s]*(x_G1[s]**2+y_G1[s]**2)\
             +I_G2[s]+M_G2[s]*(x_G2[s]**2+y_G2[s]**2)


    ##G0,G1,G2に働く力
    N_G1[s] = M_G1[s]*g*sin(theta[s]) - M_G1[s]*y_G1[s]*w[s]**2 + M_G1[s]*x_G1[s]*dw[s]
    N_G2[s] = N_G1[s] + M_G0[s]*g*cos(theta[s]) - M_G0[s]*y_G0[s]*w[s]**2 + M_G0[s]*x_G0[s]*dw[s]

    f1[s]= c*(dx_G0[s]-dx_G1[s]) + k*(x[high_unit,s]-x_G1[s])  #上
    f2[s]= c*(dx_G0[s]-dx_G2[s]) + k*(x[low_unit,s]-x_G2[s])   #下
    if theta[s] < 0.1:
        FH[s] = f2[s]
        FN[s] = M*AllUnits*g

def impulse(unit,s):
    """
    衝撃力が働いている間のG0,G1,G2のx座標y座標及び，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s :ステップ数
    """
    print("dx_G0[s]",dx_G0[s],"x_G0[s]",x_G0[s])
    ##t_s+1における速度と座標
    dx_G0[s+1] = dx_G0[s] + (F-(f1[s]+f2[s]))*dt/M
    x_G0[s+1]  = x_G0[s] + dx_G0[s]*dt
    print("dx_G0[s+1]",dx_G0[s+1],"x_G0[s+1]",x_G0[s+1])
    ##上 G1###
    ##G1の速度，座標

    dx_G1[s+1] = dx_G1[s] +f1[s]*dt/M_G1[s]
    x_G1[s+1]  = x_G1[s] + dx_G1[s]*dt

    ##角速度，角度
    dw[s] = (F*h*(2*unit+1) + M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s]))/I_G[s]
    w[s+1] = w[s]+ dw[s]*dt
    theta[s+1] =  theta[s] + w[s]*dt

    dx_G2[s+1] = 0
    x_G2[s+1] = -b
    # dw[s] = 0
    # w[s+1] = 0
    # theta[s+1] = 0

def PosRotate(unit,s):
    """
    重心座標，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s:ステップ数
    """
    ##G0の加速度
    DDXG0 = (f1[s] + f2[s])/M_G0[s] - g*sin(theta[s]) -x_G0[s]*w[s]**2 + y_G0[s]*dw[s]
    ##G0の速度
    dx_G0[s+1] = dx_G0[s] + DDXG0*dt
    ##G0のx座標
    x_G0[s+1]  = x_G0[s]  + dx_G0[s]*dt

    ##G1の加速度
    DDXG1 =  - f1[s]/M_G1[s] - g*sin(theta[s]) -x_G1[s]*w[s]**2 - y_G1[s]*dw[s]
    ##G1の速度
    dx_G1[s+1] = dx_G1[s] + DDXG1*dt
    ##G1のx'座標
    x_G1[s+1]  = x_G1[s]  + dx_G1[s]*dt

    ##角速度，角度
    dw[s]  = M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s])/I_G[s]
    w[s+1] = w[s] + dw[s]*dt
    theta[s+1] =  theta[s] + w[s]*dt

def impulse_high(s):
    """
    運動量保存則
        M_G0[s]*dx_G0[s] + M*dx_G1[s] = (M_G0[s]+M)*dx_G0[s+1]

    引数
        s:ステップ数
    """
    dx_G0[s+1] = (M_G0[s]*dx_G0[s]+M*dx_G1[s])/(M_G0[s]+M)
    x_G0[s+1]  = x_G0[s] + dx_G0[s]*dt

    DDXG1 =  - f1[s]/M_G1[s] - g*sin(theta[s]) -x_G1[s]*w[s]**2 - y_G1[s]*dw[s]
    ##G1の速度
    dx_G1[s+1] = dx_G1[s] + DDXG1*dt
    ##G1のx'座標
    x_G1[s+1]  = x_G1[s]  + dx_G1[s]*dt

    x_G2[s+1] = x_G2[s]
    ##角速度，角度
    dw[s]  = M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s])/I_G[s]
    w[s+1] = w[s] + dw[s]*dt
    theta[s+1] =  theta[s] + w[s]*dt


def impulse_low(s):
    """
    運動量保存則
        M_G0[s]*dx_G0[s] + M*dx_G2[s] = (M_G0[s]+M)*dx_G0[s+1]
    dx_G2[s] = 0より
        M_G0[s]*dx_G0[s] = (M_G0[s]+M)*dx_G0[s+1]

    引数
        s:ステップ数
    """
    dx_G0[s+1] = dx_G0[s]*M_G0[s]/(M_G0[s]+M)
    x_G0[s+1]  = x_G0[s] + dx_G0[s]*dt

    DDXG1 =  - f1[s]/M_G1[s] - g*sin(theta[s]) -x_G1[s]*w[s]**2 - y_G1[s]*dw[s]
    ##G1の速度
    dx_G1[s+1] = dx_G1[s] + DDXG1*dt
    ##G1のx'座標
    x_G1[s+1]  = x_G1[s]  + dx_G1[s]*dt


    x_G2[s+1] = x_G2[s]

    ##角速度，角度
    dw[s]  = M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s])/I_G[s]
    w[s+1] = w[s] + dw[s]*dt
    theta[s+1] =  theta[s] + w[s]*dt

def xy(low_unit,high_unit,s):

    x[:low_unit,s] = -b
    x[low_unit:high_unit+1,s] = x[low_unit:high_unit+1,s-1] + x_G0[s] - x_G0[s-1]
    x[high_unit+1:,s] = x_G1[s]
    y[:,s] = y[:,0]

    xf[:,s] = xf[:,s-1] + x[:,s]*cos(theta[s]) - y[:,s]*sin(theta[s])
    yf[:,s] = yf[:,s-1] + x[:,s]*sin(theta[s]) + y[:,s]*cos(theta[s])

low = Attacked_unit
high  = Attacked_unit
print("#####################################################################################")
print("#####################################################################################")
print("#####################################################################################")


for t in range(step-1):
    print("---------------------------------------------------------------------------------")
    print("t:",t,"low:",low,"high:",high)
    if t == 0:
        Inertia(low,high,t)
        print(0.00)
    else:
        IF1F2(low,high,t)
        print(0.1)

    if t < int(TF/dt):
        impulse(Attacked_unit,t)
        xy(low,high,t)
        print("<1>")
    else:
        if x[low,t]-(-b) > x_max:
            impulse_low(t)
            xy(low,high,t)
            low -= 1
            Inertia(low,high,t)
            print("<2>")
        elif x[high,t]-x_G1[t] > x_max:
            impulse_high(t)
            xy(low,high,t)
            high  += 1
            Inertia(low,high,t)
            print("<3>")
        else:
            PosRotate(Attacked_unit,t)
            xy(low,high,t)
            print("<4>")

    if high == AllUnits or low == 1:
        break
    if theta[t] < 0:
        theta[t] = 0

    print("x_G",x_G[t],"y_G",y_G[t],"I_G",I_G[t])
    print("x[low+1,t]",x[low+1,t],"x[high,t]",x[high,t])
    print("M_G0",M_G0[t],"x_G0[t]",x_G0[t],"dx_G0[t]",dx_G0[t])
    print("M_G1",M_G1[t],"x_G1[t]",x_G1[t],"dx_G1[t]",dx_G1[t])
    print("theta",theta[t])
    print("f1",f1[t],"f2",f2[t])
#
plt.figure(figsize=(8,3))
# t vs x のグラフ
plt.plot(x_G0,label = "xG0")
plt.plot(x_G1,label = "xG1")
plt.plot(x_G2,label = "xG2")
# #
# plt.plot(I_G0,label = "IG0")
# plt.plot(I_G1,label = "IG1")
# plt.plot(I_G2,label = "IG2")
#
# plt.plot(f1,label = "f1")
# plt.plot(f2,label = "f2")

plt.xlim(0,100)
plt.xlabel('t (step)')
plt.ylabel('x_G')
plt.legend()
plt.show()
