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
TF = 1.0     #[s] 衝撃力がかかっている時間
dt = 0.01    #[s] 微小時間
step = int(T/dt) #ステップ数


##各変数の定義
##np.array(下から数えたユニットの番号，ステップ数)
x = np.zeros( ((AllUnits,step)))    ##ユニットiのx'座標: x[i][j]，
y = np.zeros ((AllUnits,step))    ##y'座標y[i][j]

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
y[:,0] += [i*2*h+h for i in range(AllUnits)]

##[N]衝撃力
F = 10
Attacked_unit = 5

##簡略化のための変数．初期化
THETA = theta[0]
W = w[0]
DW = dw[0]
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
NG1  = N_G1[0]
NG2  = N_G2[0]
hG0  = h_G0[0]
hG1  = h_G1[0]
hG2  = h_G2[0]
F1   = f1[0]
F2   = f2[0]
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
    y_G0[s] = (y[high_unit+1,s]-y[low_unit+1,s])/2#y座標
    M_G0[s] = M*num
    h_G0[s] = h*num##2h*num/2
    I_G0[s] = 0
    for i in range(low_unit, high_unit+1):
        I_G0[s] += I + M*((x[i,s]-x_G0[s])**2+(y[i,s]-y_G0[s])**2)

    ##G1の重心の座標，質量
    x_G1[s] = sum(x[high_unit+1:, s])/(AllUnits-high_unit-1) #x座標
    y_G1[s] = sum(y[high_unit+1:, s])/(AllUnits-high_unit-1)#y座標
    M_G1[s] = M*(AllUnits-high_unit)
    h_G1[s] = M*(AllUnits-high_unit)
    I_G1[s] = M_G1[s]*((AllUnits-high_unit)**2+b**2)/3
    ##G2の重心の座標，質量
    x_G2[s] = sum(x[:low_unit, s])/low_unit #x座標
    y_G2[s] = sum(y[:low_unit, s])/low_unit#y座標
    M_G2[s] = M*low_unit
    h_G2[s] = h*low_unit
    I_G2[s] = M_G2[s]*(low_unit**2+b**2)/3


    ##G0,G1,G2に働く力
    N_G1[s] = M_G1[s]*g*sin(theta[s]) - M_G1[s]*y_G1[s]*w[s]**2 + M_G1[s]*x_G1[s]*dw[s]
    N_G2[s] = N_G1[s] + M_G0[s]*g*cos(theta[s]) - M_G0[s]*y_G0[s]*w[s]**2 + M_G0[s]*x_G0[s]*dw[s]

    f1[s]= c*(dx_G0[s]-dx_G1[s]) + k*(x[high_unit,s]-x_G1[s])  #上
    f2[s]= c*(dx_G0[s]-dx_G2[s]) + k*(x[low_unit,s]-x_G2[s])   #下
    if theta[s] < 0.1:
        FH[s] = f2[s]
        FN[s] = N_G2[s]


    #置き換え
    THETA = theta[s]
    W = w[s]
    DW = dw[s]
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

    ##下 G2##
    ##G2の速度，座標##
    dx_G2[s+1] = dx_G2[s+1] + F2*dt/M_G2[s]
    x_G2[s+1]  = XG2 + DXG2*dt
    x[:unit,s+1]=x_G2[s+1] ##G2全てのユニットのx座標をxG2に

def PosRotate(unit,s):
    """
    重心座標，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s:ステップ数
    """
    ##G0の加速度
    DDXG0 = (F1 + F2)/M_G0[s] - g*sin(THETA) -XG0*W**2 -YG0*DW
    ##G0の速度
    dx_G0[s+1] = DXG0 + DDXG0*dt
    ##G0のx座標
    x_G0[s+1]  = XG0  + DXG0*dt

    ##G1の加速度
    DDXG1 =  - F1/M_G1[s] - g*sin(THETA) -XG1*W**2 - YG1*DW
    ##G1の速度
    dx_G0[s+1] = DXG1 + DDXG1*dt
    ##G1のx'座標
    x_G0[s+1]  = XG1  + DXG1*dt

    ##N,H,dWを求めるための行列
    A = np.zeros((3,3))
    B = np.zeros(3)
    NHDW = np.zeros(3)
    ##NHDW = [N,H,dW]
    A = np.array([
    [cos(THETA), sin(THETA), M_G2[s]*YG2],
    [sin(THETA), -cos(THETA), M_G2[s]*b],
    [-YG2,-b,-IG0+IG1+IG2]])
    B = np.array([
        [M_G2[s]*(g*sin(THETA)-b*W**2)],
        [M_G2[s]*(-g*cos(THETA)+YG2*W**2)-NG2],
        [NG1*(XG0-XG1) - NG2*(XG0-XG2) + F1*(hG0-hd) + F2*(hG0+hd)]
    ])
    ##NHDW = A^(-1).b
    NHDW = np.dot(np.linalg.inv(A),B)
    FN[s] = NHDW[0]
    FH[s] = NHDW[1]
    dw[s+1] = NHDW[2]
    w[s+1] = W + DW*dt
    theta[s+1] = THETA * W*dt

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
            impulse_low()
            low -= 1
        elif x[high,t]-x_G1[t] > x_max:
            impulse_high()
            high  += 1
        else:
            PosRotate(Attacked_unit,t)
    xy(low,high,t)

plt.figure(figsize=(8,3))
# t vs x のグラフ
plt.plot(x_G0)
plt.xlabel('t (step)')
plt.ylabel('x (mm)')
plt.legend()
plt.show()




# def Inertia_low(low_unit,high_unit,s):
#     """
# 	G0,G1,G2の重心座標，質量，慣性モーメント，垂直抗力，を求める関数
#     引数
#         low_unit：フリーユニット（下）
#         high_unit：フリーユニット（上）
#         s：ステップ数
#         G0，G1，G2の重心，質量，慣性モーメントを求める
#     """
#     #間のユニットの数
#     num = high_unit - low_unit + 1
#
#     ##G0の重心の座標，質量
#     x_G0[s] = sum(x[low_unit:high_unit+1, s])/num #x座標
#     y_G0[s] = (y[high_unit+1,s]-y[low_unit+1,s])/2#y座標
#     M_G0[s] = M*num
# 	h_G0[s] = h*num	 ##2h*num/2
#
# 	##更新
# 	I_G0[s] = 0
#     for i in range(low_unit, high_unit+1):
#         I_G0[s] += I + M*((x[i,s]-x_G0[s])**2+(y[i,s]-y_G0[s])**2)
#
#     ##G1の重心の座標，質量
#     x_G1[s] = sum(x[high_unit+1:, s])/(AllUnits-high_unit-1) #x座標
#     y_G1[s] = sum(y[high_unit+1:, s])/(AllUnits-high_unit-1)#y座標
#     M_G1[s] = M*(AllUnits-high_unit)
# 	h_G1[s] = M*(AllUnits-high_unit)
#     I_G1[s] = M_G1[s]*((AllUnits-high_unit)**2+b**2)/3
#     ##G2の重心の座標，質量
#     x_G2[s] = sum(x[:low_unit, s])/low_unit #x座標
#     y_G2[s] = sum(y[:low_unit, s])/low_unit#y座標
#     M_G2[s] = M*low_unit
# 	h_G2[s] = h*low_unit
#     I_G2[s] = M_G2[s]*(low_unit**2+b**2)/3
#
#     ##G0,G1,G2に働く力
#     N_G1[s] = M_G1[s]*g*sin(theta[s]) - M_G1[s]*y_G1[s]*w[s]**2 + M_G1[s]*x_G1[s]*dw[units,s]
#     N_G2[s] = N_G1[s] + M_G0[s]*g*cos(theta[s]) - M_G0[s]*y_G0[s]*w[s]**2 + M_G0[s]*x_G0[s]*dw[s]
#
#     f1[s]= c*(dx[high_unit,s]-dx_G1[s]) + k*(x[high_unit,s]-x_G1[s])  #上
#     f2[s]= c*(dx[low_unit,s]-dx_G2[s]) + k*(x[low_unit,s]-x_G2[s])   #下
#     if theta[s] < 0.1:
#         H[s] = f2[s]
#         N[s] = N_G2[s]
#
#
#     #置き換え
#     X = x[unit,s]
#     Y = y[unit,s]
#     DX = dx[unit,s]
# 	DDX = ddx[unit,s]
#     THETA = theta[s]
#     W = w[s]
#     DW = dw[s]
# 	XG0 = x_G0[s]
#     YG0 = y_G0[s]
#     XG1 = x_G1[s]
#     YG1 = y_G1[s]
#     XG2 = x_G2[s]
#     YG2 = y_G2[s]
#     DXG0 = dx_G0[s]
#     DXG1 = dx_G1[s]
#     DXG2 = dx_G2[s]
#     DDXG0 = ddx_G0[s]
#     DDXG1 = ddx_G1[s]
#     DDXG2 = ddx_G2[s]
#     NG1  = N_G1[unit,s]
#     NG2  = N_G2[unit,s]
#     nG1  = n_G1[unit,s]
#     nG2  = n_G2[unit,s]
# 	hG0  = h_G0[s]
# 	hG1  = h_G1[s]
# 	hG2  = h_G2[s]
#     F1   = f1[s]
#     F2   = f2[s]

# def impulse_low(low_unit,s):
#     """
#     衝撃力をFとして先のFを更新する
#     F  =
#     fb =
#     """
#     ##G0について
#     dx_G0[s+1] = DXG0 + dt*(-F - (f1+f2)/M_G0[s] + sin(THETA)+XG0*W**2-YG0*DW)
#     x_G0[s+1]  = XG0  + DXG0*dt
#     ##G1について
#     dx_G1[s+1] = DXG1 + dt*(f1/M_G1[s] + g*sin(THETA)+XG1*W**2+YG1*DW)
#     x_G1[s+1]  = XG1  + DXG1*dt
#     ##A(low_unit)について
#     dx[low_unit,s+1] = dx[low_unit,s] + dt*((F+f-fb)_M_G2[s] + sin(THETA)+x[low_unit,s]*W**2 - y[low_unit,s]*DW)
#     x[low_unit,s+1]  = x[low_unit,s]  + dx[low_unit,s]*dt
#
#     fb = c*(dx[low_unit,s]-dx_G2[s]) + k*(x[last_unit,s]-x_G2[s])  #上
#
#     ##N,H,dWを求めるための行列
#     A = np.zeros(3,3)
#     b = np.zeros(3)
#     NHDW = np.zeros(3)
#     ##NHDW = [N,H,dW]
#     A = np.array([
#         [cos(THETA), sin(THETA), M_G2[s]*YG2],
#         [sin(THETA), cos(THETA), M_G2[s]*b],
#         [-YG2, -b, -IG0+IG1+IG2]
#     ])
#     b = np.array([
#         [M_G2[s]*(g*sin(THETA)-b*W**2)],
#         [M_G2[s]*(-g*cos(TEHTA)+YG2*w**2)-NG2],
#         [NG1*(XG0-XG1) - NG2*(XG0-XG2) + F1*(hG0-hG2) + F2*(hG0+hd)]
#     ])
###
