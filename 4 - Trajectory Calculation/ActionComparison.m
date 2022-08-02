clc
clear all

%%
addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\1 particle\Data')
for i = 1:19
    Name1 = ['Traj_1p_' num2str(i) '.mat'];
    data1 = load(Name1);

    data1 = data1.IterData;
    Action(i, 1) = data1.Action;
    Action(i, 2) = data1.AnalyticalAction;
    Alpha(i)        = data1.AlphaValues(i);
end

figure(1)
clf(figure(1))
hold on
plot(Alpha, Action(:, 1), 'o', 'DisplayName', 'MC action')
plot(Alpha, Action(:, 2), '.-', 'DisplayName', 'Analytica action')
xlabel('\alpha')
ylabel('S')
hold off

figure(2)
clf(figure(2))
hold on
plot(Alpha, Action(:, 1) - Action(:, 2), '.-', 'DisplayName', 'Action difference')
xlabel('\alpha')
ylabel('S diff.')
hold off

%%
clc
clear all
addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\3 particle\Standard MC Data')
for i = 1:16
    Name1 = ['Traj_3p_STDMC' num2str(i) '.mat'];
    data1 = load(Name1);

    Name2 = ['Traj_3p_RESTMC' num2str(i) '.mat'];
    data2 = load(Name2);

    data1 = data1.IterData;
    data2 = data2.IterData;

    Action(i, 1) = data1.Action;
    Action(i, 2) = data1.FittedAction;
    Action(i, 3) = data2.Action;
    Action(i, 4) = data2.FittedAction;
    Alpha(i)        = data1.AlphaValues(i);
end

figure(1)
clf(figure(1))
hold on
plot(Alpha, Action(:, 1), 'o', 'DisplayName', 'MC action')
plot(Alpha, Action(:, 3), 'ko', 'DisplayName', 'Restricted MC action')
plot(Alpha, Action(:, 2), '.-', 'DisplayName', 'Fitted MC action')
plot(Alpha, Action(:, 4), '.-', 'DisplayName', 'Restricted Fitted MC action')
xlabel('\alpha')
ylabel('S')
legend
hold off

figure(2)
clf(figure(2))
hold on
plot(Alpha, Action(:, 2) - Action(:, 4), '.-', 'DisplayName', 'Fitted Action differences')
plot(Alpha, Action(:, 1) - Action(:, 3), '.-', 'DisplayName', 'Action differences')
xlabel('\alpha')
ylabel('S diff.')
legend
hold off

%%
clc
clear all

for i = 1:16
    nameSTR = ['Traj_1p_SIMPLEMC' num2str(i)];
    Data1 = load(nameSTR);
    Data1 = Data1.IterData;

    AlphaSimple(i) = Data1.AlphaValues(i);
    ActionSimple(i) = Data1.Action;
end

for i = 1:16
    nameSTR = ['Traj_3p_RESTMC' num2str(i)];
    Data2 = load(nameSTR);
    Data2 = Data2.IterData;

    AlphaRest(i) = Data2.AlphaValues(i);
    ActionRest(i) = Data2.Action;
end

for i = 1:16
    nameSTR = ['Traj_3p_STDMC' num2str(i)];
    Data3 = load(nameSTR);
    Data3 = Data3.IterData;

    AlphaStd(i) = Data3.AlphaValues(i);
    ActionStd(i) = Data3.Action;
end

%%

figure(1)
clf(figure(1))
hold on
plot(AlphaStd, ActionStd - ActionRest, 'o')
plot(AlphaRest, ActionStd - ActionSimple, 'd')
hold off
%%
figure(2)
clf(figure(2))
hold on
plot(Data1.time, Data1.Trajectories(2, :) - Data2.Trajectories(2, :))
plot(Data2.time, Data1.Trajectories(2, :)  - Data3.Trajectories(2, :))
hold off