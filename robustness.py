#!/usr/bin/env python3
from math import *
import numpy as np
from matplotlib import pyplot as plt

##ユニットの各パラメータ
M  = 0.05 #[kg] ユニットの質量
c  = 0.1    #粘性定数
k  = 0.1    #バネ定数
b  = 200.0/2*0.001  #[m] ユニットの幅/2
h  = 30.0/2 *0.001   #[m] ユニットの高さ/2
hd = 20.0  *0.001        #[m] バネダンパの作用点までの距離
I  = M*(b*b + h*h)/3  #１ユニットの重心まわりの慣性モーメント
g  = 9.8    ##[m/s2]重力加速度
x_max = 30.0*0.001    ##最大変位量
AllUnits = 20         #ユニット数

##[N]衝撃力
F = 100
Attacked_unit = 10

##時間
T  = 5    #[s] 全体の時間
TF = 0.01     #[s] 衝撃力がかかっている時間
dt = 0.001    #[s] 微小時間
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
n_N  = np.zeros (step)  ##垂直効力のx座標
FN   = np.zeros(step)
FH   = np.zeros(step)

##初期値
x[:,0] = -b
x_G0[0] = -b
x_G1[0] = -b
x_G2[:] = -b
y[:,0] += [i*2*h+h for i in range(AllUnits)]

n_N[0] = -b

##終了時間
EndTime = 0


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
    ##x_G2 = -b
    y_G2[s] = sum(y[:low_unit, s])/low_unit#y座標
    M_G2[s] = M*low_unit
    h_G2[s] = h*low_unit
    I_G2[s] = M_G2[s]*(h_G2[s]**2+b**2)/3

    ##全体の重心位置
    x_G[s] = (M_G0[s]*x_G0[s] + M_G1[s]*x_G1[s] + M_G2[s]*x_G2[s])/(M*AllUnits)
    y_G[s] = h*AllUnits
    ##左下まわりの全体の慣性モーメント
    I_G[s] = I_G0[s]+M_G0[s]*(x_G0[s]**2+y_G0[s]**2)\
             +I_G1[s]+M_G1[s]*(x_G1[s]**2+y_G1[s]**2)\
             +I_G2[s]+M_G2[s]*(x_G2[s]**2+y_G2[s]**2)

    ##G0,G1,G2に働く力
    N_G1[s] = M_G1[s]*g*cos(theta[s]) - M_G1[s]*y_G1[s]*w[s]**2 + M_G1[s]*x_G1[s]*dw[s]
    N_G2[s] = N_G1[s] + M_G0[s]*g*cos(theta[s]) - M_G0[s]*y_G0[s]*w[s]**2 + M_G0[s]*x_G0[s]*dw[s]

    f1[s]= c*(dx_G0[s]-dx_G1[s]) + k*(x[high_unit,s]-x_G1[s])  #上
    f2[s]= c*(dx_G0[s]-dx_G2[s]) + k*(x[low_unit,s]-x_G2[s])   #下

    FH[s] = F
    FN[s] = M*AllUnits*g

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

    # print("##dx_G1[s]",dx_G1[s],"x_G1[s]",x_G1[s],"y_G1[s]",y_G1[s],"M_G1[s]",M_G1[s],"I_G1[s]",I_G1[s])
    # print("##dx_G0[s]",dx_G0[s],"x_G0[s]",x_G0[s],"y_G0[s]",y_G0[s],"M_G0[s]",M_G0[s],"I_G0[s]",I_G0[s])
    # print("##dx_G2[s]",dx_G2[s],"x_G2[s]",x_G2[s],"y_G2[s]",y_G2[s],"M_G2[s]",M_G2[s],"I_G2[s]",I_G2[s])

    ##全体の重心位置
    x_G[s] = (M_G0[s]*x_G0[s]+M_G1[s]*x_G1[s]+M_G2[s]*x_G2[s])/(M*AllUnits)
    y_G[s] = h*AllUnits

    f1[s]= c*(dx_G0[s]-dx_G1[s]) + k*(x[high_unit,s]-x_G1[s])  #上
    f2[s]= c*(dx_G0[s]-dx_G2[s]) + k*(x[low_unit,s]-x_G2[s])   #下

    FH[s] = f2[s]
    FN[s] = M*AllUnits*g

def impulse(unit,s):
    """
    衝撃力が働いている間のG0,G1,G2のx座標y座標及び，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s :ステップ数
    """

    ##t_s+1における速度と座標
    dx_G0[s+1] = dx_G0[s] + (F-(f1[s]+f2[s]))*dt/M
    x_G0[s+1]  = x_G0[s] + dx_G0[s]*dt
    ##上 G1###
    ##G1の速度，座標
    dx_G1[s+1] = dx_G1[s] +f1[s]*dt/M_G1[s]
    x_G1[s+1]  = x_G1[s] + dx_G1[s]*dt

    ##角度が0以上のとき，角速度，角度を計算する．
    n_N[s] = -(F*h*(2*unit+1)+M*AllUnits*g*x_G[s])/FN[s]
    if x_G[s] > 0 or theta[s] > 0 or n_N[s] > -0.01:
        ##角速度，角度
        n_N[s] = 0  ##垂直抗力のx座標は0
        dw[s] = (F*h*(2*unit+1) + M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s]))/I_G[s]
        w[s+1] = w[s] + dw[s]*dt
        theta[s+1] =  theta[s] + w[s]*dt
    else:
        n_N[s] = x_G[s]

def PosRotate(unit,s):
    """
    重心座標，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s:ステップ数
    """
    ##G0の加速度
    DDXG0 = -(f1[s] + f2[s])/M_G0[s] + g*sin(theta[s])  + x_G0[s]*w[s]**2 - y_G0[s]*dw[s]
    ##G0の速度
    dx_G0[s+1] = dx_G0[s] + DDXG0*dt
    ##G0のx座標
    x_G0[s+1]  = x_G0[s]  + dx_G0[s]*dt

    ##G1の加速度
    DDXG1 =   f1[s]/M_G1[s] + g*sin(theta[s]) + x_G1[s]*w[s]**2 - y_G1[s]*dw[s]
    ##G1の速度
    dx_G1[s+1] = dx_G1[s] + DDXG1*dt
    ##G1のx'座標
    x_G1[s+1]  = x_G1[s]  + dx_G1[s]*dt

    ##角度が0以上のとき，角速度，角度を計算する．
    if x_G[s] > 0 or theta[s] > 0:
        ##角速度，角度
        n_N[s] = 0  ##垂直抗力のx座標は0
        dw[s]  = M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s])/I_G[s]
        w[s+1] = w[s] + dw[s]*dt
        theta[s+1] =  theta[s] + w[s]*dt
    else:
        n_N[s] = x_G[s]

def impulse_high(low_unit,high_unit,unit,s):
    """
    運動量保存則
        M_G0[s]*dx_G0[s] + M*dx_G1[s] = (M_G0[s]+M)*dx_G0[s+1]

    引数
        s:ステップ数
    """
    dx_G0[s+1] = (M_G0[s]*dx_G0[s]+M*dx_G1[s])/(M_G0[s]+M)
    Inertia(low_unit,high_unit+1,s)
    x_G0[s+1]  = x_G0[s] + dx_G0[s]*dt

    DDXG1 =  - f1[s]/(M_G1[s]-M) + g*sin(theta[s]) + x_G1[s]*w[s]**2 - y_G1[s]*dw[s]
    ##G1の速度
    dx_G1[s+1] = dx_G1[s] + DDXG1*dt
    ##G1のx'座標
    x_G1[s+1]  = x_G1[s]  + dx_G1[s]*dt

    ##角度が0以上のとき，角速度，角度を計算する．
    if x_G[s] > 0 or theta[s] > 0:
        ##角速度，角度
        n_N[s] = 0  ##垂直抗力のx座標は0
        dw[s]  = M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s])/I_G[s]
        w[s+1] = w[s] + dw[s]*dt
        theta[s+1] =  theta[s] + w[s]*dt
    else:
        n_N[s] = x_G[s]

def impulse_low(low_unit,high_unit,unit,s):
    """
    運動量保存則
        M_G0[s]*dx_G0[s] + M*dx_G2[s] = (M_G0[s]+M)*dx_G0[s+1]
    dx_G2[s] = 0より
        M_G0[s]*dx_G0[s] = (M_G0[s]+M)*dx_G0[s+1]

    引数
        s:ステップ数
    """
    dx_G0[s+1] = dx_G0[s]*M_G0[s]/(M_G0[s]+M)
    Inertia(low_unit-1,high_unit,s)
    x_G0[s+1]  = x_G0[s] + dx_G0[s]*dt

    DDXG1 =  - f1[s]/M_G1[s] + g*sin(theta[s]) + x_G1[s]*w[s]**2 - y_G1[s]*dw[s]
    ##G1の速度
    dx_G1[s+1] = dx_G1[s] + DDXG1*dt
    ##G1のx'座標
    x_G1[s+1]  = x_G1[s]  + dx_G1[s]*dt


    ##角度が0以上のとき，角速度，角度を計算する．
    if x_G[s] > 0 or theta[s] > 0:
        ##角速度，角度
        n_N[s] = 0  ##垂直抗力のx座標は0
        dw[s]  = M*AllUnits*g*(sin(theta[s])*y_G[s] + cos(theta[s])*x_G[s])/I_G[s]
        w[s+1] = w[s] + dw[s]*dt
        theta[s+1] =  theta[s] + w[s]*dt
    else:
        n_N[s] = x_G[s]


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

flag = 1

for t in range(step-1):
    #print("---------------------------------------------------------------------------------")

    if t == 0:
        Inertia(low,high,t)

    else:
        IF1F2(low,high,t)

    ##衝撃
    if t < int(TF/dt):
        impulse(Attacked_unit,t)
        xy(low,high,t)
        print("<1>")
    ##衝撃後
    else:
        if x[low,t-1]-(-b) > x_max:
            impulse_low(low,high,Attacked_unit,t)
            low -= 1
            xy(low,high,t)
            print("t:",t*dt,"low:",low,"high:",high)

            print("<2>")

        elif x[high,t-1]-x_G1[t-1] > x_max:
            impulse_high(low,high,Attacked_unit,t)
            high  += 1
            xy(low,high,t)
            print("t:",t*dt,"low:",low,"high:",high)

            print("<3>")

        else:
            PosRotate(Attacked_unit,t)
            xy(low,high,t)
            #print("<4>")

    if high == 18 or low == 1:
        print(EndTime)
        break

    if theta[t] > pi/2:
        theta[t+1:] = theta[t]
        print(EndTime)
        break

    EndTime = t

#
plt.figure(figsize=(8,5))
# t vs x のグラフ

for i in range (AllUnits):
    plt.plot(x[i],label = i)
plt.xlabel('t [ms]')
plt.ylabel('x[m]')
plt.legend()
plt.xlim(0,EndTime)
plt.savefig('各ユニットのx座標.png')
plt.show()

plt.plot(x_G,label = "xG")
plt.plot(x_G0,label = "xG0")
plt.plot(x_G1,label = "xG1")
plt.plot(x_G2,label = "xG2")
plt.xlabel('t [ms]')
plt.ylabel('x[m]')
plt.legend()
plt.xlim(0,EndTime)
plt.savefig('全体G0G1G2の重心.png')
plt.show()

plt.plot(dx_G0,label = "dxG0")
plt.plot(dx_G1,label = "dxG1")
plt.plot(dx_G2,label = "dxG2")
plt.xlabel('t [ms]')
plt.ylabel('dxdt[m/s]')
plt.legend()
plt.xlim(0,EndTime)
plt.savefig('重心速度.png')
plt.show()


plt.plot(f1,label = "f1")
plt.plot(f2,label = "f2")
plt.xlabel('t [ms]')
plt.ylabel('f1')
plt.legend()
plt.xlim(0,EndTime)
plt.show()

plt.plot(theta,label = "theta")
# plt.plot(w,label = "w")
plt.xlabel('t [ms]')
plt.ylabel('theta')
plt.legend()
plt.xlim(0,EndTime)
plt.savefig('角度_吸収あり.png')
plt.show()

#
# plt.plot(I_G,label = "IG")
# plt.xlabel('t [ms]')
# plt.ylabel('IG')
# plt.legend()
# plt.xlim(0,EndTime)
# plt.show()

plt.plot(n_N,label = "n")
plt.xlabel('t [ms]')
plt.ylabel('n')
plt.legend()
plt.xlim(0,EndTime)
plt.savefig('垂直抗力の位置.png')
plt.show()


#######################################################################
#######################################################################

dw2 = np.zeros(step)
w2 = np.zeros(step)
theta2 = np.zeros(step)
n2 = np.zeros(step)

M2 = M*AllUnits
xG2 = -b
yG2 = h*AllUnits

I2 = 4*M2*(yG2**2+b**2)/3
H2 = h*(2*Attacked_unit+1)
n2[:] = xG2
EndTime2 = 0
for t in range(step-1):
    if t < int(TF/dt):
        n2[t] = (F*H2 + M2*g*(sin(theta2[t])*y_G[t] + cos(theta2[t])*x_G[t]))/(M2*g)
        print(t,xG2,n2[t],theta2[t],F*H2,M2*g*xG2)
        if n2[t] > 0:
            ##角速度，角度
            n2[t] = 0  ##垂直抗力のx座標は0
            dw2[t]  = (F*H2 + M2*g*(sin(theta2[t])*yG2 + cos(theta2[t])*xG2))/I2
            w2[t+1] = w2[t] + dw2[t]*dt
            theta2[t+1] =  theta2[t] + w2[t]*dt

    else:
        print(t,xG2,n2[t],theta2[t])
        if theta2[t] > 0:
            ##角速度，角度
            n2[t] = 0  ##垂直抗力のx座標は0
            dw2[t]  = M2*g*(sin(theta2[t])*yG2 + cos(theta2[t])*xG2)/I2
            w2[t+1] = w2[t] + dw2[t]*dt
            theta2[t+1] =  theta2[t] + w2[t]*dt
        else:
            n2[t] = xG2

        if theta2[t] == 0:

            break

    if theta2[t] > pi/2:
        theta2[t+1:] = theta2[t]
        break

    EndTime2 = t

plt.figure(figsize=(8,5))
plt.plot(theta,label = "theta")
plt.plot(theta2,label = "theta2")
# plt.plot(w,label = "w")
plt.xlabel('t [ms]')
plt.ylabel('theta')
plt.legend()
if EndTime > EndTime2:
    plt.xlim(0,EndTime)
else:
    plt.xlim(0,EndTime2)
plt.savefig('角度.png')
plt.show()
#
# plt.figure(figsize=(8,3))
# # t vs x のグラフ
# plt.plot(theta2,label = "theta2")
#
# # #
# # plt.plot(I_G0,label = "IG0")
# # plt.plot(I_G1,label = "IG1")
# # plt.plot(I_G2,label = "IG2")
# #
# # plt.plot(f1,label = "f1")
# # plt.plot(f2,label = "f2")
#
# plt.xlim(0,1000)
# plt.xlabel('t (step)')
# plt.ylabel('theta2')
# plt.legend()
# plt.show()
