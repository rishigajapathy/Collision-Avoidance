% Solving ODE using Runge-Kutta 4th order method
clear 
load XsimCopy.mat;   
h=1;
t_os = 0:h:3500;
N = length(t_os);

u_os = zeros(1,N);
v_os = zeros(1,N); 
r_os = zeros(1,N);
psi_os = zeros(1,N);
x_os = zeros(1,N);
y_os = zeros(1,N);


Init_cond_os = [7.97 0 0 0 0 0];
p_rps_os = 1.783;
rud_def_i_os = 35*pi/180;

u_os(1) = Init_cond_os(1);
v_os(1) = Init_cond_os(2);
r_os(1) = Init_cond_os(3);
psi_os(1) = Init_cond_os(4);
x_os(1) = Init_cond_os(5);
y_os(1) = Init_cond_os(6);


rdot = zeros(1,N);
rdot(1) = 0;



for i=1:N-1

    [due1,dve1,dre1,dpsie1,dxe1,dye1]= dyn_os(t_os(i),u_os(i),v_os(i),r_os(i),psi_os(i),rdot(i),p_rps_os,rud_def_i_os);

    rdot(i+1)= dre1;
    
    k1 = h*due1; l1 = h*dve1; m1 = h*dre1; n1 = h*dpsie1; o1 = h*dxe1; p1 = h*dye1;
                        
    [due2,dve2,dre2,dpsie2,dxe2,dye2]= dyn_os(t_os(i)+0.5*h,u_os(i)+0.5*k1,v_os(i)+0.5*l1,r_os(i)+0.5*m1,psi_os(i)+0.5*n1,rdot(i),p_rps_os,rud_def_i_os);

    k2 = h*due2; l2 = h*dve2; m2 = h*dre2; n2 = h*dpsie2; o2 = h*dxe2; p2 = h*dye2;
    
    [due3,dve3,dre3,dpsie3,dxe3,dye3]= dyn_os(t_os(i)+0.5*h,u_os(i)+0.5*k2,v_os(i)+0.5*l2,r_os(i)+0.5*m2,psi_os(i)+0.5*n2,rdot(i),p_rps_os,rud_def_i_os);
 
    k3 = h*due3; l3 = h*dve3; m3 = h*dre3; n3 = h*dpsie3; o3 = h*dxe3; p3 = h*dye3;

    [due4,dve4,dre4,dpsie4,dxe4,dye4]= dyn_os(t_os(i)+h,u_os(i)+k3,v_os(i)+l3,r_os(i)+m3,psi_os(i)+n3,rdot(i),p_rps_os,rud_def_i_os);

    k4 = h*due4; l4 = h*dve4; m4 = h*dre4; n4 = h*dpsie4; o4 = h*dxe4; p4 = h*dye4;

    u_os(i+1) = u_os(i) + (k1+2*k2+2*k3+k4)/6; 
    v_os(i+1) = v_os(i) + (l1+2*l2+2*l3+l4)/6;
    r_os(i+1) = r_os(i) + (m1+2*m2+2*m3+m4)/6;
    psi_os(i+1) = psi_os(i) + (n1+2*n2+2*n3+n4)/6;
    x_os(i+1) = x_os(i) + (o1+2*o2+2*o3+o4)/6;
    y_os(i+1) = y_os(i) + (p1+2*p2+2*p3+p4)/6;
 
end


plot(XsimCopy(:,5),XsimCopy(:,4))
hold on
plot(y_os,x_os)
hold off
axis equal
grid on







