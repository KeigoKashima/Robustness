#!/usr/bin/env python3
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
dx = np.zeros ((AllUnits,step))   ##xの変化量dx/dt
ddx = np.zeros ((AllUnits,step))  ##二階微分 ddx/dt
theta = np.zeros ((AllUnits,step))##角度
w = np.zeros ((AllUnits,step))    ##角速度
dw = np.zeros ((AllUnits,step))   ##角加速度

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


def Inertia(first_unit,last_unit,s):
    """
	G0,G1,G2の重心座標，質量，慣性モーメント，垂直抗力，を求める関数
    引数
        first_unit：フリーユニット（下）
        last_unit：フリーユニット（上）
        s：ステップ数
        G0，G1，G2の重心，質量，慣性モーメントを求める
    """
    #間のユニットの数
    num = last_unit - first_unit + 1

    ##G0の重心の座標，質量
    x_G0[s] = sum(x[first_unit:last_unit+1, s])/num #x座標
    y_G0[s] = (y[last_unit+1,s]-y[first_unit+1,s])/2#y座標
    M_G0[s] = M*num
	h_G0[s] = h*num	 ##2h*num/2

	##更新
	I_G0[s] = 0
    for i in range(first_unit, last_unit+1):
        I_G0[s] += I + M*((x[i,s]-x_G0[s])**2+(y[i,s]-y_G0[s])**2)

    ##G1の重心の座標，質量
    x_G1[s] = sum(x[last_unit+1:, s])/(AllUnits-last_unit-1) #x座標
    y_G1[s] = sum(y[last_unit+1:, s])/(AllUnits-last_unit-1)#y座標
    M_G1[s] = M*(AllUnits-last_unit)
	h_G1[s] = M*(AllUnits-last_unit)
    I_G1[s] = M_G1[s]*((AllUnits-last_unit)**2+b**2)/3
    ##G2の重心の座標，質量
    x_G2[s] = sum(x[:first_unit, s])/first_unit #x座標
    y_G2[s] = sum(y[:first_unit, s])/first_unit#y座標
    M_G2[s] = M*first_unit
	h_G0[s] = h*first_unit
    I_G2[s] = M_G2[s]*(first_unit**2+b**2)/3

    ##G0,G1,G2に働く力
    N_G1[s] = M_G1[s]*g*sin(theta[s]) - M_G1[s]*y_G1[s]*w[unit,s]**2 + M_G1[s]*x_G1[s]*dw[units,s]
    N_G2[s] = N_G1[s] + M_G0[s]*g*cos(theta[s]) - M_G0[s]*y_G0[s]*w[unit,s]**2 + M_G0[s]*x_G0[s]*dw[unit,s]
    N[s]    =

def impulse(unit,s):
    """
	衝撃力が働いている間のG0,G1,G2のx座標y座標及び，回転角度を求める関数
    引数
        unit:衝撃が与えられたユニット
        s :ステップ数
    """

    #置き換え
    X = x[unit,s]
    Y = y[unit,s]
    DX = dx[unit,s]
    THETA = theta[unit,s]
    W = w[unit,s]
    DW = dw[unit,s]
    XG1 = x_G1[s]
    YG1 = y_G1[s]
    XG2 = x_G2[s]
    YG2 = y_G2[s]
    DXG1 = dx_G1[s]
    DXG2 = dx_G2[s]
    #DYG1 = dy_G1[s]
    #DYG2 = dy_G2[s]
    NG1  = N_G1[unit,s]
    NG2  = N_G2[unit,s]
    nG1  = n_G1[unit,s]
    nG2  = n_G2[unit,s]


    ##衝撃力がかかっている間の各変数の値を求める
    ##ユニット間 G0###
    f1[unit,s]= c*(DX-DXG1) + k*(X-XG1)  #上
    f2[unit,s]= c*(DX-DXG2) + k*(X-XG2)  #下
    ##t_s+1における速度と座標
    dx[unit,s+1] = DX + (F-(f1[unit,s]+f2[unit,s]))*dt/M
    x[unit,s+1]  = X + DX*dt
    ##t_s+1における角速度と角度
    w[unit,s+1] = W + ((NG1*nG1+NG2*nG2) + f1[unit,s]*(h-hd) + f2[unit,s]*(h+hd))*dt/IG0
    theta[unit,s+1]=THETA+W*dt

    ##上 G1###
    ##G1の速度，座標
    dx_G1[s+1] = dx_G1[s+1] + f1[unit,s]*dt/M_G1[s]
    x_G1[s+1]  = XG1 + DXG1*dt
    x[unit+1:,s+1]=x_G1[s+1] ##G1全てのユニットのx座標をxG1に
    ##G1の角速度角度
    #w[unit,s+1]= (-f1[unit,s]*(y_G1[s]-Y-hd) + NG1*(nG1-X+XG1))*dt/I_G1
    #theta[unit,s+1] = THETA+W*dt

    ##下 G2##
    ##G2の速度，座標##
    dx_G2[s+1] = dx_G2[s+1] + f2[unit,s]*dt/M_G2[s]
    x_G2[s+1]  = XG2 + DXG2*dt
    x[:unit,s+1]=x_G2[s+1] ##G2全てのユニットのx座標をxG2に
    ##G2の角速度角度
    #w[unit,s+1]= (f2[unit,s]*(Y-y_G2[s]+hd)+NG2*(nG2+X-XG2))*dt/I_G2
    #theta[unit,s+1] = THETA+W*dt

def G0(s):
    """
	G0の重心座標，回転角度を求める関数
    引数
        s:ステップ数
    """
    #置き換え
    X = x[unit,s]
    Y = y[unit,s]
    DX = dx[unit,s]
	DDX = ddx[unit,s]
    THETA = theta[unit,s]
    W = w[unit,s]
    DW = dw[unit,s]
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
    DYG1 = dy_G1[s]
    DYG2 = dy_G2[s]
    DDXG0 = ddx_G0[s]
    DDXG1 = ddx_G1[s]
    DDXG2 = ddx_G2[s]
    DDYG0 = ddy_G0[s]
    DDYG1 = ddy_G1[s]
    DDYG2 = ddy_G2[s]
    NG1  = N_G1[unit,s]
    NG2  = N_G2[unit,s]
    nG1  = n_G1[unit,s]
    nG2  = n_G2[unit,s]
	hG0  = h_G0[s]
	hG1  = h_G1[s]
	hG2  = h_G2[s]
	f1   = f1[unit,s]
	f2   = f2[unit,s]


    ##G0の加速度
    ddx_G0[s+1] = DDXG0 + (c*(2*DXG0-dx[last+1,s]-dx[first-1,s]) + k*(x[last,s]-x[last+1,s] + x[first,s]-x[first-1,s]))/M_G0[s] \
					- g*sin(THETA) -XG0*W**2 -YG0*DW
    ##G0の速度
    dx_G0[unit,s+1] = DXG1 + DDXG1*dt
    ##G0のx座標
    x_G0[unit,s+1]  = XG1  + DXG1*dt

    ##G0の角加速度
	dw[unit,s] = DW + (NG1*nG1 + NG2*nG2 +  f1*(hG0-hd) + f2*(hG0+hd))*dt
    ##G0の角速度
    w_G0[unit,s+1] = WG0 + DWG1*dt
    ##G0の角度
    theta_G0[unit,s+1]  = THETA + W*dt

def G1(s):
    """
	G1の重心座標，回転角度を求める関数
    引数
        s:ステップ数
    """
    #置き換え
    X = x[unit,s]
    Y = y[unit,s]
    DX = dx[unit,s]
	DDX = ddx[unit,s]
    THETA = theta[unit,s]
    W = w[unit,s]
    DW = dw[unit,s]
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
    DYG1 = dy_G1[s]
    DYG2 = dy_G2[s]
    DDXG0 = ddx_G0[s]
    DDXG1 = ddx_G1[s]
    DDXG2 = ddx_G2[s]
    DDYG0 = ddy_G0[s]
    DDYG1 = ddy_G1[s]
    DDYG2 = ddy_G2[s]
    NG1  = N_G1[unit,s]
    NG2  = N_G2[unit,s]
    nG1  = n_G1[unit,s]
    nG2  = n_G2[unit,s]
	hG0  = h_G0[s]
	hG1  = h_G1[s]
	hG2  = h_G2[s]
	f1   = f1[unit,s]
	f2   = f2[unit,s]

    ##G1の加速度
    ddx_G1[unit,s+1] = DDXG1 - f1/M_G1[s] - g*sin(THETA) -XG1*W**2 - YG1*DW
    ##G1の速度
    dx_G0[unit,s+1] = DXG1 + DDXG1*dt
    ##G1のx'座標
    x_G0[unit,s+1]  = XG1  + DXG1*dt

    ##G1の角加速度
	dw[unit,s] = DW + NG1*(nG1-XG0+XG1) - f1*(hG1+hd)
    ##G1の角速度
	w_G1[unit,s+1] = WG1 + DWG1*dt
    ##G1の角度
    theta_G1[unit,s+1]  = THETA + W*dt

def G2(s):
    """
	G2の重心座標，回転角度を求める関数
    引数
        s:ステップ数
    """
    #置き換え
    X = x[unit,s]
    Y = y[unit,s]
    DX = dx[unit,s]
	DDX = ddx[unit,s]
    THETA = theta[unit,s]
    W = w[unit,s]
    DW = dw[unit,s]
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
    DYG1 = dy_G1[s]
    DYG2 = dy_G2[s]
    DDXG0 = ddx_G0[s]
    DDXG1 = ddx_G1[s]
    DDXG2 = ddx_G2[s]
    DDYG0 = ddy_G0[s]
    DDYG1 = ddy_G1[s]
    DDYG2 = ddy_G2[s]
    NG1  = N_G1[unit,s]
    NG2  = N_G2[unit,s]
    nG1  = n_G1[unit,s]
    nG2  = n_G2[unit,s]
	hG0  = h_G0[s]
	hG1  = h_G1[s]
	hG2  = h_G2[s]
	f1   = f1[unit,s]
	f2   = f2[unit,s]

    ##G2の加速度
    ddx_G2[unit,s+1] = DDXG2 - (f2 - N*sin(THETA) - H*cos(THETA))/M_G2[s] \
                       - g*sin(THETA) - XG2*W**2 + YG2*DW
    ##G2の速度
    dx_G2[unit,s+1] = DXG2 + DDXG2*dt
    ##G2のx'座標
    x_G2[unit,s+1]  = XG2  + DXG2*dt

    ##G2の角加速度
	dw[unit,s] = DW + NG2*(nG2+XG0-XG2) + N*XG2 + H*YG2 - f2*(hG2-hd)
    ##G2の角速度
	w_G2[unit,s+1] = WG2 + DWG2*dt
    ##G2の角度
    theta_G2[unit,s+1]  = THETA + W*dt



first = Attacked_unit
last  = Attacked_unit

for t in range(int(TF/dt)):
    Inertia(first,last,t)

    if t<TF/dt:
        impulse(Attacked_unit,t)

    else:
        if x[first,t]-x_G2[t] > x_max:
            first -= 1
        elif x[last,t]-x_G1[t] > x_max:
            last  += 1
        else:
            G0(t)
        #    G1(t)
        #    G2(t)

print(x_G0[0:10])
print(x_G1[0:10])
print(x_G2[0:10])


plt.figure(figsize=(8,3))
# t vs x のグラフ
plt.plot(x_G0)
plt.xlim(0, 10)
plt.ylim(-50,-49.5)
plt.xlabel('t (step)')
plt.ylabel('x (m)')
plt.show()
