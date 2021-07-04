
% ==============  Scientific Computing - Project 2  =================

%%
approximation_u10 = importdata("approximation_u10.txt");
approximation_u100 = importdata("approximation_u100.txt");
approximation_u200 = importdata("approximation_u200.txt");
% Plotting
figure(1)
plot(approximation_u10(:,2),approximation_u10(:,1),"LineWidth",1.0)
hold on 
plot(approximation_u100(:,2),approximation_u100(:,1),"-.","LineWidth",1.0,"color","r")
hold on 
plot(approximation_u200(:,2),approximation_u200(:,1),"--","LineWidth",1.0,"color","k")
grid on
legend("N=10","N=100","N=200")
xlabel("x")
ylabel("u")
%% 
% Graph of solution to advection equation
x = linspace(0,2*pi,251);
C=1;
t = [0,0.25,0.5,1];
figure(2)
for i=1:length(t)
    t(i);
    u = 1.5 + sin(x-C*t(i));
    plot(x,u)
    hold on
end
xlabel("x")
ylabel("u")
xlim([0 7])
ylim([0 3])
legend("t=0","t=0.25","t=0.5","t=1","Location","Best")
grid on
% Solution obtained from Galerkin method
time0 = importdata("time0.txt");
time025 = importdata("time025.txt");
time05 = importdata("time05.txt");
time1 = importdata("time1.txt");
x = linspace(0,2*pi,length(time0));
figure(3)
plot(x,time0,"-.","LineWidth",0.9)
hold on
plot(x,time025,"--","LineWidth",0.9)
hold on
plot(x,time05,":","color","k","LineWidth",0.9)
hold on
plot(x,time1,"LineWidth",0.9)
xlabel("x")
ylabel("u")
legend("t=0","t=0.25","t=0.5","t=1")
grid on 
%% Grid and timestep resolution
x = linspace(0,2*pi,101);
C=1;
t = [0,0.25,0.5,1];
figure(8)
for i=1:length(t)
    t(i);
    u = 1.5 + sin(x-C*t(i));
    p1 = plot(x,u,"color","b");
    hold on
end
time0 = importdata("time0.txt");
time025 = importdata("time025.txt");
time05 = importdata("time05.txt");
time1 = importdata("time1.txt");
x = linspace(0,2*pi,length(time0));
plot(x,time0,"-.","color","k")
hold on
plot(x,time025,"-.","color","k")
hold on
plot(x,time05,"-.","color","k")
hold on
p2 = plot(x,time1,"-.","color","k");
xlabel("x")
ylabel("u")
legend([p1 p2],{"Analytic solution","Numerical solution"})
%% Question 2
uQ2_0 = importdata("Q2data0.txt");
uQ2 = importdata("Q2data.txt");
x = linspace(0,2*pi,length(uQ2_0));
uQ2_1 = importdata("Q2data1.txt");
uQ2_2 = importdata("Q2data2.txt");
figure(1)
plot(x,uQ2)
grid on
xlabel("x")
ylabel("u")
ylim([-0.2 1.2])
%
figure(2)
plot(x,uQ2_0)
grid on
xlabel("x")
ylabel("u")
%
figure(3)
plot(x,uQ2_1)
grid on
xlabel("x")
ylabel("u")
%
figure(4)
plot(x,uQ2_2)
grid on
xlabel("x")
ylabel("u")
%legend("t=0","t=0.25","t=0.5","t=1")
% Highly oscillatory at the discontinuities, not at all physical, or what
% is supposed to happen with a discontinuous intial velocity profile.

%% Question 3
% Solving Burger's equation using 
u0 = importdata("Q3data0.txt");
u7 = importdata("Q3data7.txt");
u9 = importdata("Q3data9.txt");
u11 = importdata("Q3data11.txt");
x = linspace(0,2*pi,length(u7));
figure(5)
plot(x,u7,"-.")  
hold on
plot(x,u9,"-.")
hold on
plot(x,u0)
hold on
plot(x,u11,"-.")
grid on 
xlabel("x")
ylabel("u")
%% Sine-wave IC 
uQ3_0 = importdata("Q3sine_0.txt");
uQ3_05 = importdata("Q3sine_05.txt");
uQ3_09 = importdata("Q3sine_09.txt");
uQ3_1 = importdata("Q3sine_1.txt");
uQ3_11 = importdata("Q3sine_11.txt");
uQ3_125 = importdata("Q3sine_125.txt");
uQ3_15 = importdata("Q3sine_15.txt");
uQ3_175 = importdata("Q3sine_175.txt");
uQ3_2 = importdata("Q3sine_2.txt");
x = linspace(0,2*pi,length(uQ3_0));
figure()
plot(x,uQ3_0)
hold on
plot(x,uQ3_05)
hold on
plot(x,uQ3_09)
ylim([0 3])
xlabel("x")
ylabel("u")
grid on

figure()
plot(x,uQ3_1)
hold on
plot(x,uQ3_11)
xlabel("x")
ylabel("u")
grid on


figure()
plot(x,uQ3_125)
hold on
plot(x,uQ3_15)
xlabel("x")
ylabel("u")
grid on


figure()
plot(x,uQ3_175)
hold on
plot(x,uQ3_2)
hold on
xlabel("x")
ylabel("u")
grid on
%% Square-wave initial condition
uQ3_0 = importdata("Q3square0.txt");
uQ3_025 = importdata("Q3square025.txt");
uQ3_075 = importdata("Q3square075.txt");
uQ3_1 = importdata("Q3square1.txt");
uQ3_125 = importdata("Q3square125.txt");
uQ3_175= importdata("Q3square175.txt");
uQ3_2= importdata("Q3square2.txt");
x = linspace(0,2*pi,length(uQ3_2));

figure()
plot(x,uQ3_0)
xlabel("x")
ylabel("u")
ylim([-0.2 1.2])
grid on
hold on
plot(x,uQ3_025)
xlabel("x")
ylabel("u")
grid on

figure()
plot(x,uQ3_075)
xlabel("x")
ylabel("u")

hold on
plot(x,uQ3_1)
xlabel("x")
ylabel("u")
grid on

figure()
plot(x,uQ3_125)
xlabel("x")
ylabel("u")
grid on

hold on
plot(x,uQ3_175)
xlabel("x")
ylabel("u")
grid on

figure()
plot(x,uQ3_2)
xlabel("x")
ylabel("u")
grid on
