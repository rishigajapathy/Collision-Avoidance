function [dudt,dvdt,drdt,dpsidt,dxdt,dydt] = dyn_os(~,u,v,r,psi,rdot,p_rps_os,rud_def_i_os)

u_os = u;v_os = v;r_os = r; psi_os = psi;

  
%%%%%%%%%%%%% same as dynamics function code start %%%%%%%%%%%%%%%%%%%%%%%%

% mass dimesionaless mx' my' iz'
mxp=0.022;
myp=0.223;
jzp=0.011;

% characteristic dimensions - Fullscale
Lpp = 320;
d = 20.8;
% b = 58;
% cb = 0.81;
rho=1025;

% mass normalization
ndm=0.5*rho*d*Lpp^2;
ndj=0.5*rho*d*Lpp^4;

% A matrix terms 
m=312600*rho;
mx=mxp*ndm;
my=myp*ndm;
xg=11.2;
izg=1.9900e+12;
jz=jzp*ndj;

% to avoid NaN
eps=0;%1e-16;
U=sqrt(u_os^2+v_os^2+eps);
% nondimensional states u' v' r' and their powers
vp=v_os/(U+eps);
rp=r_os*Lpp/(U+eps);
vp2=v_os^2/(U^2+eps);
vp3=v_os^3/(U^3+eps);
vp4=v_os^4/(U^4+eps);
rp2=r_os^2*Lpp^2/(U^2+eps);
rp3=r_os^3*Lpp^3/(U^3+eps);

% hydrodynamic coefficients nondimensional
Xh=-0.022-0.04*vp2+0.002*vp*rp+0.011*rp2+0.771*vp4;
Yh=-0.315*vp+0.083*rp-1.607*vp3+0.379*vp2*rp-0.391*vp*rp2+0.008*rp3;
Nh=-0.137*vp-0.049*rp-0.03*vp3-0.294*vp2*rp+0.055*vp*rp2-0.013*rp3;

% hydrodynamic force
Fh=0.5*rho*U^2*Lpp*d*[Xh
    Yh
    Nh*Lpp];

% control input np-propellor speed & del-rudder angle
np=p_rps_os;
del = rud_def_i_os;


% thrust model
k2=-0.1385;
k1=-0.2753;
k0=0.2931;
Dp=9.86;
wp0=0.35;
xpp=-0.48;
bet=atan(-v_os/(u_os+eps));
betp=bet-xpp*rp;
tp=0.220;

wp=wp0*exp(-4*betp^2);
Jp=(u_os*(1-wp))/(np*Dp+eps);
KT=k0+k1*Jp+k2*Jp^2;
T=(1-tp)*rho*(np^2)*(Dp^4)*(KT);
% Thrust force
Ft=[T;0;0];

% rudder model

tr=0.387;
ah=0.312;
xhp=-0.464;
%xh=xhp*Lpp;
xr=-0.5*Lpp;
fa=2.747;
Dp=9.86;
Hr=15.8;
eta=Dp/Hr;
epsi=1.09;
kappa=0.5;
lrp=-0.710;
Ar=112.5;

ur=epsi*u_os*(1-wp)*sqrt( eta*( 1+ kappa*( sqrt(1+8*KT/(pi*Jp^2+eps))-1) )^2+(1-eta) );
betr=bet-lrp*rp;

a=0.395;
b=0.64;
gamr=(b-a)/2*tanh(betr*50)+(a+b)/2;
vr=U*gamr*betr;
alpr=del-atan(vr/(ur+eps));
Ur=sqrt(ur^2+vr^2+eps);
Fn=0.5*rho*Ar*Ur^2*fa*sin(alpr);
% rudder force
Fr=[-(1-tr)*Fn*sin(del)
    -(1+ah)*Fn*cos(del)
    -(xr+ah*xhp)*Fn*cos(del)];

%%%%%%%%%%%%% dynamics function code end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X = Fh(1) + Ft(1) + Fr(1);
Y = Fh(2) + Ft(2) + Fr(2);
Nm = Fh(3) + Ft(3) + Fr(3);
 


%%% the following are the differential equations being solved in ode45 
dudt =  ((m+my)*v_os*r_os + xg*m*r_os^2 + X)/(m+mx);
dvdt = ((-m-mx)*u_os*r_os -(m*xg*rdot)+Y)/(m+my);
drdt =  (Nm-(xg*m*(dvdt+u_os*r_os)))/(izg+xg^2*m+jz);
dpsidt = r_os;
dxdt =  u_os*cos(psi_os)-v_os*sin(psi_os);                        
dydt =  u_os*sin(psi_os)+v_os*cos(psi_os);                        
% dxdt =  u_os*sin(psi_os)+v_os*cos(psi_os);                        
% dydt =  u_os*cos(psi_os)-v_os*sin(psi_os);


end