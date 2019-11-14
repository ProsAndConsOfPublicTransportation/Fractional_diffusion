aM=0.9999175;
delta0=2;
deltaK=0.5;
mash_max_size=0.01;
initial_u=0;
t=linspace(0,1000,1000);
bc_fun=@(region,state)500;

k=800;
rho = 8960;
specificHeat = 386;
thick = .01;
stefanBoltz = 5.670373e-8;
hCoeff = 0;
ta = 0;
emiss = 0;
width=1.5;
height=1.5;
ww=width/6;
hh=height/6;
x1=ww; x2=ww*3; x3=ww*5; y1=hh*1; y2=hh*3; y3=hh*5; 
gam=1+(delta0*width/18)^2;
c = sprintf('(-(floor(0.5*floor(%g-(x-0.75).^2)+0.5*floor(%g-(y-0.75).^2))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2))))*%g+1)*%g*%g', gam,gam,gam,x1,gam, y1,gam, x1,gam, y2,gam, x1,gam, y3,gam, x2,gam, y1,gam, x2,gam, y3,gam, x3,gam, y1,gam, x3,gam, y2,gam, x3,gam, y3,aM, thick, k);
a = sprintf('2*%g + 2*%g*%g*u.^3', hCoeff, emiss, stefanBoltz);
f = 2*hCoeff*ta + 2*emiss*stefanBoltz*ta^4;
d = thick*rho*specificHeat;

gdmTrans = [3 4 0  width width 0 0 0 height height;
           3 4 3*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3;
           3 4 ww-ww*delta0/3  ww+ww*delta0/3 ww+ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww+ww*delta0/3 ww+ww*delta0/3;
           3 4 ww-ww*delta0/3  ww+ww*delta0/3 ww+ww*delta0/3 ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3;
           3 4 ww-ww*delta0/3  ww+ww*delta0/3 ww+ww*delta0/3 ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww+ww*delta0/3 5*ww+ww*delta0/3;
           3 4 5*ww-ww*delta0/3  5*ww+ww*delta0/3 5*ww+ww*delta0/3 5*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3;
           3 4 5*ww-ww*delta0/3  5*ww+ww*delta0/3 5*ww+ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww+ww*delta0/3 5*ww+ww*delta0/3;
           3 4 5*ww-ww*delta0/3  5*ww+ww*delta0/3 5*ww+ww*delta0/3 5*ww-ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww+ww*delta0/3 ww+ww*delta0/3;
           3 4 3*ww-ww*delta0/3  3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww+ww*delta0/3 5*ww+ww*delta0/3;
           3 4 3*ww-ww*delta0/3  3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww+ww*delta0/3 ww+ww*delta0/3;];

sf2=['R1+R2+R3+R4+R5+R6+R7+R8+R9+P1'];
ns2=[82 82 82 82 82 82 82 82 82 80; 49 50 51 52 53 54 55 56 57 49];
geometry = decsg(gdmTrans',sf2,ns2);

[u1,pb]=FEM_plate_sim(c,a,f,d,geometry,initial_u,t,bc_fun,mash_max_size,true,0,false);
u1_upper_edge=Temperature_at_edge(u1,1.5,0.75,pb);

figure;
    pdeplot(pb, 'xydata', u1(:,end)/10, 'contour', 'on', 'colormap', 'jet');

%FIND_FIRST_ORDER
[alpha_a_1,error_1]=Find_order(specificHeat,k,rho,u1_upper_edge);

%SIMULATE_SYSTEM_AFTER_SWITCH
delta0=deltaK;
gam=1+(delta0*width/18)^2;
c = sprintf('(-(floor(0.5*floor(%g-(x-0.75).^2)+0.5*floor(%g-(y-0.75).^2))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2)))+floor(0.5*floor(abs(%g-(x-%g).^2))+0.5*floor(abs(%g-(y-%g).^2))))*%g+1)*%g*%g', gam,gam,gam,x1,gam, y1,gam, x1,gam, y2,gam, x1,gam, y3,gam, x2,gam, y1,gam, x2,gam, y3,gam, x3,gam, y1,gam, x3,gam, y2,gam, x3,gam, y3,aM, thick, k);
a = sprintf('2*%g + 2*%g*%g*u.^3', hCoeff, emiss, stefanBoltz);
f = 2*hCoeff*ta + 2*emiss*stefanBoltz*ta^4;
d = thick*rho*specificHeat;

gdmTrans = [3 4 0  width width 0 0 0 height height;
           3 4 3*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3;
           3 4 ww-ww*delta0/3  ww+ww*delta0/3 ww+ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww+ww*delta0/3 ww+ww*delta0/3;
           3 4 ww-ww*delta0/3  ww+ww*delta0/3 ww+ww*delta0/3 ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3;
           3 4 ww-ww*delta0/3  ww+ww*delta0/3 ww+ww*delta0/3 ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww+ww*delta0/3 5*ww+ww*delta0/3;
           3 4 5*ww-ww*delta0/3  5*ww+ww*delta0/3 5*ww+ww*delta0/3 5*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww-ww*delta0/3 3*ww+ww*delta0/3 3*ww+ww*delta0/3;
           3 4 5*ww-ww*delta0/3  5*ww+ww*delta0/3 5*ww+ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww+ww*delta0/3 5*ww+ww*delta0/3;
           3 4 5*ww-ww*delta0/3  5*ww+ww*delta0/3 5*ww+ww*delta0/3 5*ww-ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww+ww*delta0/3 ww+ww*delta0/3;
           3 4 3*ww-ww*delta0/3  3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww-ww*delta0/3 5*ww+ww*delta0/3 5*ww+ww*delta0/3;
           3 4 3*ww-ww*delta0/3  3*ww+ww*delta0/3 3*ww+ww*delta0/3 3*ww-ww*delta0/3 ww-ww*delta0/3 ww-ww*delta0/3 ww+ww*delta0/3 ww+ww*delta0/3;];

sf2=['R1+R2+R3+R4+R5+R6+R7+R8+R9+P1'];
ns2=[82 82 82 82 82 82 82 82 82 80; 49 50 51 52 53 54 55 56 57 49];
geometry = decsg(gdmTrans',sf2,ns2);
initial_u=u1;
t=linspace(1000,2000,1000);

[u2,pb2]=FEM_plate_sim(c,a,f,d,geometry,initial_u,t,bc_fun,mash_max_size,false,pb,false);
u2_upper_edge=Temperature_at_edge(u2,1.5,0.75,pb2);

figure;
    pdeplot(pb2, 'xydata', u2(:,end)/10, 'contour', 'on', 'colormap', 'jet');
    
%JOINED_SOLUTIONS
u_all=[u1_upper_edge,u2_upper_edge];

figure;
    plot(u_all)

%CHECK_SECOND_ORDER
initial_u=0;
[u2_0,pb2_0]=FEM_plate_sim(c,a,f,d,geometry,initial_u,t,bc_fun,mash_max_size,true,0,false);
u2_upper_edge_0=Temperature_at_edge(u2_0,1.5,0.75,pb2_0);
[alpha_a_2,error_2]=Find_order(specificHeat,k,rho,u2_upper_edge_0);



a0=alpha_a_1(2);
a1=alpha_a_2(2);
alpha0=alpha_a_1(1);
alpha1=alpha_a_2(1);