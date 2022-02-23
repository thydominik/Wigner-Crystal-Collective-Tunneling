clc
clear all

%This script will contain all calculations needed to get the Tunneling
%energy splitting. (3 particle only)

%Constants for the calculations ------------------------------------------
PartNum     = 3;                %Number of Particles in the calculation
AlphaValNum = 30;               %Number of Alpha values 

AlphaMin    = -6;               %Minimum value of abs(alpha)
AlphaMax    = -15;              %Maximum value of abs(alpha)

Alpha       = linspace(AlphaMin, AlphaMax, AlphaValNum);

eta         = 20;               %interaction strenght (40: Pascu Calculations, 20: Upper bound for my Calculations, 18.8...:The fitted value)

EqPosIter   = 1 * 10^5;         %Number of MC iterations for the equilibrium position
TrajsIter   = 5 * 10^6;             %Number of MC iterations for the trajectories

TrajDiv     = 100;              %How many points ina a trajectory
Z           = linspace(-1, 1, TrajDiv);     %Z time -> reparametrized tau
dz          = 2 * Z(TrajDiv/2 + 1);         %time division


%Constants for the calculations ------------------------------------------

%Part 1: A MC calculation for the equilibrium positions:
[EquilibriumPositions, DiscardedEqPos, FinalEnergy] = EqPosMC(PartNum, Alpha, eta, EqPosIter);
%%
figure(2)
clf(figure(2))
hold on
plot(Alpha, EquilibriumPositions', '.-')

hold off
%%
%Part 2: MC for the trajectories using part 1:
[Trajectories, DiscardedTrajs, FinalEnergyTrajs, R_param] = TrajsMC(EquilibriumPositions, eta, Alpha, TrajsIter, Z, dz, TrajDiv);

%%
%Part 3: Prefactor Calculations
Results = PrefactorCalculation(PartNum, TrajDiv, R_param, Alpha, eta, Z, Trajectories)


