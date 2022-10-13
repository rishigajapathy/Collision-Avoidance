import numpy as np
import matplotlib.pyplot as plt

import math 

######### dyn_os start ###############

def dyn_os(u_os,v_os,r_os,psi_os,rdot_os,p_rps_os,rud_action_os):
    
    
#mass dimesionaless mx' my' iz'
  mxp=0.022
  myp=0.223
  jzp=0.011

#characteristic dimensions - Fullscale
  Lpp = 320
  d = 20.8

  rho=1025

#mass normalization
  ndm=0.5*rho*d*Lpp**2
  ndj=0.5*rho*d*Lpp**4

# A matrix terms 
  m=312600*rho
  mx=mxp*ndm
  my=myp*ndm
  xg=11.2
  izg=1.9900*10**12
  jz=jzp*ndj
  # to avoid NaN
  eps=0;#1e-16;
  U = np.sqrt((u_os**2)+(v_os**2))
#nondimensional states u' v' r' and their powers
  vp=v_os/(U+eps)
  rp=r_os*Lpp/(U+eps)
  vp2=v_os**2/(U**2+eps)
  vp3=v_os**3/(U**3+eps)
  vp4=v_os**4/(U**4+eps)
  rp2=r_os**2*Lpp**2/(U**2+eps)
  rp3=r_os**3*Lpp**3/(U**3+eps)
  # hydrodynamic coefficients nondimensional
  Xh=-0.022-0.04*vp2+0.002*vp*rp+0.011*rp2+0.771*vp4
  Yh=-0.315*vp+0.083*rp-1.607*vp3+0.379*vp2*rp-0.391*vp*rp2+0.008*rp3
  Nh=-0.137*vp-0.049*rp-0.03*vp3-0.294*vp2*rp+0.055*vp*rp2-0.013*rp3
    
    
  # hydrodynamic force
  Fh=([[0.5*rho*U**2*Lpp*d*Xh],[0.5*rho*U**2*Lpp*d*Yh],[0.5*rho*U**2*Lpp*d*Nh*Lpp]])
  

#control input np-propellor speed & del-rudder angle
  nprop=p_rps_os
  delta = rud_action_os   
  #thrust model
  k2=-0.1385
  k1=-0.2753
  k0=0.2931
  Dp=9.86
  wp0=0.35
  xpp=-0.48
  bet=math.atan(-v_os/(u_os+eps))
  betp=bet-xpp*rp
  tp=0.220

  wp=wp0*math.exp(-4*betp**2)
  Jp=(u_os*(1-wp))/(nprop*Dp+eps)
  KT=k0+k1*Jp+k2*Jp**2
  T=(1-tp)*rho*(nprop**2)*(Dp**4)*(KT)
# Thrust force
  Ft=np.array([[T],[0],[0]])
  # rudder model

  tr=0.387
  ah=0.312
  xhp=-0.464
  #xh=xhp*Lpp
  xr=-0.5*Lpp
  fa=2.747
  Dp=9.86
  Hr=15.8
  eta=Dp/Hr
  epsi=1.09
  kappa=0.5
  lrp=-0.710
  Ar=112.5
 
  ur=epsi*u_os*(1-wp)*np.sqrt(eta*( 1+ kappa*(np.sqrt(1+8*KT/(math.pi*Jp**2+eps))-1) )**2+(1-eta))
  betr=bet-lrp*rp  
  a=0.395
  b=0.64
  gamr=(b-a)/2*math.tanh(betr*50)+(a+b)/2
  vr=U*gamr*betr
  alpr=delta-np.arctan(vr/(ur+eps))
  Ur=np.sqrt(ur**2+vr**2+eps);
  Fn=0.5*rho*Ar*Ur**2*fa*np.sin(alpr)
    
  # rudder force
  Fr=np.array([[-(1-tr)*Fn*np.sin(delta)],[-(1+ah)*Fn*np.cos(delta)],[-(xr+ah*xhp)*Fn*np.cos(delta)]])
      
      
  
  X = Fh[0] + Ft[0] + Fr[0]
  Y = Fh[1] + Ft[1] + Fr[1]
  Nm = Fh[2] + Ft[2] + Fr[2]
    
     
    
### the following are the differential equations being solved in rk method 
  dudt =  ((m+my)*v_os*r_os + xg*m*r_os**2 + X)/(m+mx)
  dvdt = ((-m-mx)*u_os*r_os -(m*xg*rdot_os)+Y)/(m+my);
  drdt =  (Nm-(xg*m*(dvdt+u_os*r_os)))/(izg+xg**2*m+jz);
  dpsidt = r_os;
  dxdt =  u_os*np.cos(np.rad2deg(psi_os))-v_os*np.sin(np.rad2deg(psi_os));                        
  dydt =  u_os*np.sin(np.rad2deg(psi_os))+v_os*np.cos(np.rad2deg(psi_os));

  d_by_dt = [[dudt],[dvdt],[drdt],[dpsidt],[dxdt],[dydt]] # 6x1 
    
    
    
    
  return d_by_dt

######### dyn_os end ###############

h=1
time = 3500
time_span = np.arange(0,time+1)
N = len(time_span)

u_os = np.zeros((1,N))
v_os = np.zeros((1,N))
r_os = np.zeros((1,N))
psi_os = np.zeros((1,N))
x_os = np.zeros((1,N))
y_os = np.zeros((1,N))

#Init_cond_os = array([7.97,0,0,0,0,0])
Init_cond_os = np.array([[7.97],[0],[0],[np.deg2rad(0)],[0],[0]])
Init_cond_os_dash = np.transpose(Init_cond_os)
p_rps_os = 1.783
rud_def_i_os = np.deg2rad(35)

u_os[0,0] = Init_cond_os_dash[0,0]
v_os[0,0] = Init_cond_os_dash[0,1]
r_os[0,0] = Init_cond_os_dash[0,2]
psi_os[0,0] = Init_cond_os_dash[0,3]
x_os[0,0] = Init_cond_os_dash[0,4]
y_os[0,0] = Init_cond_os_dash[0,5]

rdots = list()


ss = range(0,N)

# for i=1:N-1
for i in range(0,N):
    
  u = u_os[0][i];v = v_os[0][i]; r = r_os[0][i]; psi = psi_os[0][i]; rdot = r_os[0][i]; 
  d_by_dt = dyn_os(u,v,r,psi,rdot,p_rps_os,rud_def_i_os)
  ddt = d_by_dt[2]
  rdots.append(ddt)
  k1 = h*d_by_dt[0][0];l1 = h*d_by_dt[1][0];m1 = h*d_by_dt[2][0];n1 = h*d_by_dt[3][0];
  o1 = h*d_by_dt[4][0];p1 = h*d_by_dt[5][0];
  
  u = u_os[0][i]+0.5*k1;v = v_os[0][i]+0.5*l1;r = r_os[0][i]+0.5*m1;psi = psi_os[0][i]+0.5*n1;
  d_by_dt = dyn_os(u,v,r,psi,rdot,p_rps_os,rud_def_i_os)
  k2 = h*d_by_dt[0][0]; l2 = h*d_by_dt[1][0]; m2 = h*d_by_dt[2][0]; n2 = h*d_by_dt[3][0]; 
  o2 = h*d_by_dt[4][0]; p2 = h*d_by_dt[5][0];
  
  u = u_os[0][i]+0.5*k2;v = v_os[0][i]+0.5*l2;r = r_os[0][i]+0.5*m2;psi = psi_os[0][i]+0.5*n2;
  print(v_os[0][i],"see here")
  
  d_by_dt = dyn_os(u,v,r,psi,rdot,p_rps_os,rud_def_i_os)
  k3 = h*d_by_dt[0][0]; l3 = h*d_by_dt[1][0]; m3 = h*d_by_dt[2][0]; n3 = h*d_by_dt[3][0]; 
  o3 = h*d_by_dt[4][0]; p3 = h*d_by_dt[5][0];
  
  u = u_os[0][i]+k3;v = v_os[0][i]+l3;r = r_os[0][i]+m3;psi = psi_os[0][i]+n3;
  d_by_dt = dyn_os(u,v,r,psi,rdot,p_rps_os,rud_def_i_os)
  k4 = h*d_by_dt[0][0]; l4 = h*d_by_dt[1][0]; m4 = h*d_by_dt[2][0]; n4 = h*d_by_dt[3][0]; 
  o4 = h*d_by_dt[4][0]; p4 = h*d_by_dt[5][0];
  
  

  u_os.append= u_os[-1] + (k1+2*k2+2*k3+k4)/6
  v_os.append= v_os[-1] + (l1+2*l2+2*l3+l4)/6
  r_os.append= r_os[-1] + (m1+2*m2+2*m3+m4)/6
  psi_os.append= psi_os[-1] + (n1+2*n2+2*n3+n4)/6
  x_os.append= x_os[-1] + (o1+2*o2+2*o3+o4)/6
  y_os.append= y_os[-1] + (p1+2*p2+2*p3+p4)/6
  
  plt.plot(x_os,y_os)
  
######################################


