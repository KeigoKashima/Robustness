#!/usr/bin/env python3
from math import *
import numpy as np
from matplotlib import pyplot as plt

##ユニットの各パラメータ
M  = 0.1 #[kg] ユニットの質量
c  = 1.0    #粘性定数
k  = 1.0    #バネ定数
b  = 100.0/2  #[mm] ユニットの幅/2
h  = 50.0/2    #[mm] ユニットの高さ/2
hd = 20.0          #[mm] バネダンパの作用点までの距離
I  = M*(b*b + h*h)/3  #１ユニットの重心まわりの慣性モーメント
g  = 9.8    ##[m/s2]重力加速度
x_max = 30.0    ##最大変位量
AllUnits = 10         #ユニット数

##時間
T  = 10.0    #[s] 全体の時間
TF = 0.1     #[s] 衝撃力がかかっている時間
dt = 0.01    #[s] 微小時間
step = int(T/dt) #ステップ数


##各変数の定義
##np.array(下から数えたユニットの番号，ステップ数)
x = np.zeros( ((AllUnits,step)))    ##ユニットiのx'座標: x[i][j]，
y = np.zeros ((AllUnits,step))    ##y'座標y[i][j]

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

##[N]衝撃力
F = 100
Attacked_unit = 5

##簡略化のための変数．初期化
THETA = theta[0]
W = w[0]
XG = x_G[0]
YG = y_G[0]
XG0 = x_G0[0]
YG0 = y_G0[0]
XG1 = x_G1[0]
YG1 = y_G1[0]
XG2 = x_G2[0]
YG2 = y_G2[0]
DXG0 = dx_G0[0]
DXG1 = dx_G1[0]
DXG2 = dx_G2[0]
DYG0 = dx_G0[0]
DDXG0 = ddx_G0[0]
DDXG1 = ddx_G1[0]
DDXG2 = ddx_G2[0]
NG1  = N_G1[0]
NG2  = N_G2[0]
nG1  = n_G1[0]
nG2  = n_G2[0]
hG0  = h_G0[0]
hG1  = h_G1[0]
hG2  = h_G2[0]
F1   = f1[0]
F2   = f2[0]

IG   = I_G[0]
IG0  = I_G0[0]
IG1  = I_G1[0]
IG2  = I_G2[0]



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
    M_G1[s] = M*(AllUnits-high_unit)
    h_G1[s] = h*(AllUnits-high_unit)
    I_G1[s] = M_G1[s]*((AllUnits-high_unit)**2+b**2)/3
    ##G2の重心の座標，質量
    x_G2[s] = sum(x[:low_unit, s])/low_unit #x座標
    y_G2[s] = sum(y[:low_unit, s])/low_unit#y座標
    M_G2[s] = M*low_unit
    h_G2[s] = h*low_unit
    I_G2[s] = M_G2[s]*(low_unit**2+b**2)/3

    ##全体の重心位置
    x_G[s] = (M_G0[s]*x_G0[s]+M_G1[s]*x_G1[s]+M_G2[s]*x_G2[s])/(M*AllUnits)
    y_G[s] = h*AllUnits
    ##左下まわりの全体の慣性モーメント
    I_G[s] = I_G0[s]+I_G1[s]+I_G2[s]+(x_G0[s]**2+y_G0[s]**2)+(x_G1[s]**2+y_G1[s]**2)+(x_G2[s]**2+y_G2[s]**2)


    ##G0,G1,G2に働く力
    N_G1[s] = M_G1[s]*g*sin(theta[s]) - M_G1[s]*y_G1[s]*w[s]**2 + M_G1[s]*x_G1[s]*dw[s]
    N_G2[s] = N_G1[s] + M_G0[s]*g*cos(theta[s]) - M_G0[s]*y_G0[s]*w[s]**2 + M_G0[s]*x_G0[s]*dw[s]

    f1[s]= c*(dx_G0[s]-dx_G1[s]) + k*(x[high_unit,s]-x_G1[s])  #上
    f2[s]= c*(dx_G0[s]-dx_G2[s]) + k*(x[low_unit,s]-x_G2[s])   #下
    if theta[s] < 0.1:
        FH[s] = f2[s]
        FN[s] = M*AllUnits*g


    #置き換え
    THETA = theta[s]
    W  = w[s]
    XG = x_G[s]
    YG = y_G[s]
    XG0 = x_G0[s]
    YG0 = y_G0[s]
    XG1 = x_G1[s]
    YG1 = y_G1[s]
    XG2 = x_G2[s]
    YG2 = y_G2[s]
    DXG0 = dx_G0[s]
    DXG1 = dx_G1[s]
    DXG2 = dx_G2[s]
    DYG0 = dx_G0[s]
    DDXG0 = ddx_G0[s]
    DDXG1 = ddx_G1[s]
    DDXG2 = ddx_G2[s]
    NG1  = N_G1[s]
    NG2  = N_G2[s]
    nG1  = n_G1[s]
    nG2  = n_G2[s]
    hG0  = h_G0[s]
    hG1  = h_G1[s]
    hG2  = h_G2[s]
    F1   = f1[s]
    F2   = f2[s]

    IG   = I_G[s]
    IG0  = I_G0[s]
    IG1  = I_G1[s]
    IG2  = I_G2[s]

def impulse(unit,s):
    """
    衝撃力が働いている間のG0,G1,G2のx座標y座標及び，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s :ステップ数
    """

    ##t_s+1における速度と座標
    dx_G0[s+1] = DXG0 + (F-(F1+F2))*dt/M
    x_G0[s+1]  = XG0 + DXG0*dt

    ##上 G1###
    ##G1の速度，座標
    dx_G1[s+1] = DXG1 + F1*dt/M_G1[s]
    x_G1[s+1]  = XG1 + DXG1*dt
    x[unit+1:,s+1]=x_G1[s+1] ##G1全てのユニットのx座標をxG1に

    ##角速度，角度
    dw[s] = (F*h*(2*unit+1) + M*AllUnits*g*(sin(THETA)*YG + cos(THETA)*XG))/IG
    w[s+1] = W+ dw[s]*dt
    theta[s+1] =  THETA + W*dt

def PosRotate(unit,s):
    """
    重心座標，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s:ステップ数
    """
    ##G0の加速度
    DDXG0 = (F1 + F2)/M_G0[s] - g*sin(THETA) -XG0*W**2 + YG0*dw[s]
    ##G0の速度
    dx_G0[s+1] = DXG0 + DDXG0*dt
    ##G0のx座標
    x_G0[s+1]  = XG0  + DXG0*dt

    ##G1の加速度
    DDXG1 =  - F1/M_G1[s] - g*sin(THETA) -XG1*W**2 - YG1*dw[s]
    ##G1の速度
    dx_G0[s+1] = DXG1 + DDXG1*dt
    ##G1のx'座標
    x_G0[s+1]  = XG1  + DXG1*dt

    ##角速度，角度
    dw[s]  = M*AllUnits*g*(sin(THETA)*YG + cos(THETA)*XG)/IG
    w[s+1] = W + dw[s]*dt
    theta[s+1] =  THETA + W*dt



def impulse_high(s):
    """
    運動量保存則
        M_G0[s]*DXG0 + M*DXG1 = (M_G0[s]+M)*dx_G0[s+1]

    引数
        s:ステップ数
    """
    dx_G0[s+1] = (M_G0[s]*DXG0+M*DXG1)/(M_G0[s]+M)

def impluse_low(s):
    """
    運動量保存則
        M_G0[s]*DXG0 + M*DXG2 = (M_G0[s]+M)*dx_G0[s+1]
    DXG2 = 0より
        M_G0[s]*DXG0 = (M_G0[s]+M)*dx_G0[s+1]

    引数
        s:ステップ数
    """
    dx_G0[s+1] = DXG0*M_G0[s]/(M_G0[s]+M)

def xy(low_unit,high_unit,s):
    X = np.zeros(AllUnits)
    X[:low_unit] = x_G1[s] - x_G1[s-1]
    X[low_unit:high_unit+1] = x_G0[s] - x_G0[s-1]
    X[high_unit+1:] = x_G2[s] - x_G2[s-1]
    Y = y[:,0]
    x[:,s] = x[:,s-1] + X*cos(THETA) + Y*sin(THETA)
    y[:,s] = y[:,s-1] - X*sin(THETA) + Y*cos(THETA)

low = Attacked_unit
high  = Attacked_unit

for t in range(step-1):

    Inertia(low,high,t)

    if t < int(TF/dt):
        impulse(Attacked_unit,t)

    else:
        if x[low,t]-x_G2[t] > x_max:
            impulse_low(t)
            low -= 1
        elif x[high,t]-x_G1[t] > x_max:
            impulse_high(t)
            high  += 1
        else:
            PosRotate(Attacked_unit,t)
    if t > 0:
        xy(low,high,t)

plt.figure(figsize=(8,3))
# t vs x のグラフ
plt.plot(x_G0)
plt.xlabel('t (step)')
plt.ylabel('IG')
plt.legend()
plt.show()
