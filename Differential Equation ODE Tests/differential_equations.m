clear all
clc
%1ST; 1DB LINEÁRIS, SZEPARÁLHATÓ ESET
%y(t)' = a^2 - y(t)^2

a       = 5;    %some parameter; scalar
T_st    = -7;    %time start at *here*
T_fin   = 7;   %time finishes at *here*
y0      = 0;    %some again random initial condition

[T, Y] = diff1(a,T_st,T_fin,y0);
figure(1)
clf(figure(1))
hold on
plot(T,Y,'r','LineWidth',3)
hold off

% 2ND 1DB NEMLINEÁRIS -> VÁLTOZÓ PARAMÉTERES
%y(t)' = a(t)^2 - y(t)^2
clear all
clc
a           = linspace(0, 5,1000);
T_st        = -7;
T_fin       = 7;
time        = linspace(T_st, T_fin, length(a));
y0          = 0;
[timeval, yval] = diff2(time,a,y0);

figure(2)
clf(figure(2))
hold on
plot(timeval,yval,'r','LineWidth',3)
hold off

% 3RD 2DB LINEÁRIS DIFFEGYENLET.
%y(t)' = a^2 - y(t)^2
%z(t)' = y(t) - z(t)
clear all
clc

a = 5;
T_st        = -7;
T_fin       = 7;
y0          = 0;
z0          = 0;

[T,Y] = diff3(a,T_st, T_fin, y0, z0);

figure(3)
clf(figure(3))
hold on
plot(T,Y(:,1),'r','LineWidth',3)
plot(T,Y(:,2),'k','LineWidth',3)
hold off

% LAST ONE: DIFFEGYENLET RENDSZER IDÕFÜGGÕ PARAMÉTERREL, HA EZ MEGY AKKOR MINDEN MEGY
%y(t)' = a(t)^2 - y(t)^2
%z(t)' = y(t) - z(t)
clear all
clc
a           = linspace(0, 5,1000);
T_st        = -7;
T_fin       = 7;
time        = linspace(T_st, T_fin, length(a));
y0          = 0;
z0          = 0;


[T,Y] = diff4(a, y0, z0, time);

figure(4)
clf(figure(4))
hold on
plot(T,Y(1,:),'r','LineWidth',3)
plot(T,Y(2,:),'k','LineWidth',3)
hold off









