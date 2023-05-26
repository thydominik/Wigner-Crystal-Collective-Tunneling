clear all
clc

% Constants:
Nx1 = 20;
Nx2 = 20;
Nx3 = 50;
Nx4 = Nx2;
Nx5 = Nx1;

Alpha   = 14;
Kappa   = 0;
eta = 20;

Polarization = struct();
Beta = 10^-8;

Eq_Pos = [];

Nx1 * Nx2 * Nx3 * Nx4 * Nx5

for AlphaInd = 1:length(Alpha)
    for KappaInd = 1:length(Kappa)

        % calculating the EqPos then construct the space of each particles:
        a = Alpha(AlphaInd);
        k = Kappa(KappaInd);
        Potential = @(x) 0.25 * ((x(1)^2 - a)^2 + (x(2)^2 - a)^2 + (x(3)^2 - a)^2 + (x(4)^2 - a)^2 + (x(5)^2 - a)^2) + k * x(1) + k * x(2) + k * x(3) + k * x(4) + k * x(5) + eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - x(5)) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - x(5)) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - x(5)) + 1/abs(x(4) - x(5)));
        % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
        options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
        x_start = [-sqrt(a)-1 -sqrt(a)+0.5 -sqrt(a)+1 sqrt(a)-1 sqrt(a)+1];
        [x0, fval0] = fminsearch(Potential, x_start, options);
        Eq_Pos = sort(x0);

        a = 1.5;
        XMin1   = Eq_Pos(1) - a;
        XMax1   = -Eq_Pos(5) + a;
        b = 3;
        XMin2   = Eq_Pos(2) - b;
        XMax2   = -Eq_Pos(4) + b;
        c = 3.5;
        XMin3   = Eq_Pos(3) - c;
        XMax3   = -Eq_Pos(3) + c;
        d = 3;
        XMin4   = Eq_Pos(4) - d;
        XMax4   = -Eq_Pos(2) + d;
        e = 1.5;
        XMin5   = Eq_Pos(5) - e;
        XMax5   = -Eq_Pos(1) + e;

        % Kinetic Matrix:
        x1 = linspace(XMin1, XMax1, Nx1);
        x2 = linspace(XMin2, XMax2, Nx2);
        x3 = linspace(XMin3, XMax3, Nx3);
        x4 = linspace(XMin4, XMax4, Nx4);
        x5 = linspace(XMin5, XMax5, Nx5);

        dx1 = x1(2) - x1(1);
        dx2 = x2(2) - x2(1);
        dx3 = x3(2) - x3(1);
        dx4 = x4(2) - x4(1);
        dx5 = x5(2) - x5(1);

        K1 = -1/(2 * dx1^2);
        K2 = -1/(2 * dx2^2);
        K3 = -1/(2 * dx3^2);
        K4 = -1/(2 * dx4^2);
        K5 = -1/(2 * dx5^2);

        for i = 1:Nx1
            if i == 1
                KineticMtx1(i, i)        = -2 * K1;
            else
                KineticMtx1(i - 1, i)    = 1 * K1;
                KineticMtx1(i, i - 1)    = 1 * K1;
            end
        end

        for i = 1:Nx2
            if i == 1
                KineticMtx2(i, i)        = -2 * K2;
            else
                KineticMtx2(i - 1, i)    = 1 * K2;
                KineticMtx2(i, i - 1)    = 1 * K2;
            end
        end

        for i = 1:Nx3
            if i == 1
                KineticMtx3(i, i)        = -2 * K3;
            else
                KineticMtx3(i - 1, i)    = 1 * K3;
                KineticMtx3(i, i - 1)    = 1 * K3;
            end
        end

        for i = 1:Nx4
            if i == 1
                KineticMtx4(i, i)        = -2 * K4;
            else
                KineticMtx4(i - 1, i)    = 1 * K4;
                KineticMtx4(i, i - 1)    = 1 * K4;
            end
        end

        for i = 1:Nx5
            if i == 1
                KineticMtx5(i, i)        = -2 * K5;
            else
                KineticMtx5(i - 1, i)   = 1 * K5;
                KineticMtx5(i, i - 1)   = 1 * K5;
            end
        end
        KineticMtx1 = sparse(KineticMtx1);
        KineticMtx2 = sparse(KineticMtx2);
        KineticMtx3 = sparse(KineticMtx3);
        KineticMtx4 = sparse(KineticMtx4);
        KineticMtx5 = sparse(KineticMtx5);

        K1          = kron(kron(KineticMtx1, speye(Nx2)), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5))));
        K2          = kron(speye(Nx1), kron(KineticMtx2, kron(speye(Nx3), kron(speye(Nx4), speye(Nx5)))));
        K3          = kron(speye(Nx1), kron(speye(Nx2), kron(KineticMtx3, kron(speye(Nx4), speye(Nx5)))));
        K4          = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(KineticMtx4, speye(Nx5)))));
        K5          = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(speye(Nx4), KineticMtx5))));
        KineticMtx  = K1 + K2 + K3 + K4 + K5;

        clear K1 K2 K3 K4 K5

        PotentialMtx1 = sparse(zeros(Nx1));
        PotentialMtx2 = sparse(zeros(Nx2));
        PotentialMtx3 = sparse(zeros(Nx3));
        PotentialMtx4 = sparse(zeros(Nx4));
        PotentialMtx5 = sparse(zeros(Nx5));

        U1 = 0.25 * (x1.^2 - Alpha(AlphaInd)).^2 + Kappa(KappaInd) * x1;
        U2 = 0.25 * (x2.^2 - Alpha(AlphaInd)).^2 + Kappa(KappaInd) * x2;
        U3 = 0.25 * (x3.^2 - Alpha(AlphaInd)).^2 + Kappa(KappaInd) * x3;
        U4 = 0.25 * (x4.^2 - Alpha(AlphaInd)).^2 + Kappa(KappaInd) * x4;
        U5 = 0.25 * (x5.^2 - Alpha(AlphaInd)).^2 + Kappa(KappaInd) * x5;

        PotentialMtx1 = sparse(diag(U1));
        PotentialMtx2 = sparse(diag(U2));
        PotentialMtx3 = sparse(diag(U3));
        PotentialMtx4 = sparse(diag(U4));
        PotentialMtx5 = sparse(diag(U5));

        U1 = kron(kron(PotentialMtx1, speye(Nx2)), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5))));
        U2 = kron(kron(speye(Nx1), PotentialMtx2), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5))));
        U3 = kron(kron(speye(Nx1), speye(Nx2)), kron(PotentialMtx3, kron(speye(Nx4), speye(Nx5))));
        U4 = kron(kron(speye(Nx1), speye(Nx2)), kron(speye(Nx3), kron(PotentialMtx4, speye(Nx5))));
        U5 = kron(kron(speye(Nx1), speye(Nx2)), kron(speye(Nx3), kron(speye(Nx4), PotentialMtx5)));
        PotentialMtx = U1 + U2 + U3 + U4 + U5;

        clear U1 U2 U3 U4 U5

        Hamiltonian = KineticMtx + PotentialMtx;

        X_matrix1   = sparse(diag(x1));
        X_matrix2   = sparse(diag(x2));
        X_matrix3   = sparse(diag(x3));
        X_matrix4   = sparse(diag(x4));
        X_matrix5   = sparse(diag(x5));

        X1  = kron(X_matrix1, kron(speye(Nx2), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5)))));
        X2  = kron(speye(Nx1), kron(X_matrix2, kron(speye(Nx3), kron(speye(Nx4), speye(Nx5)))));
        X3  = kron(speye(Nx1), kron(speye(Nx2), kron(X_matrix3, kron(speye(Nx4), speye(Nx5)))));
        X4  = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(X_matrix4, speye(Nx5)))));
        X5  = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(speye(Nx4), X_matrix5))));

        X12 = diag(X1) - diag(X2);
        X13 = diag(X1) - diag(X3);
        X14 = diag(X1) - diag(X4);
        X15 = diag(X1) - diag(X5);
        X23 = diag(X2) - diag(X3);
        X24 = diag(X2) - diag(X4);
        X25 = diag(X2) - diag(X5);
        X34 = diag(X3) - diag(X4);
        X35 = diag(X3) - diag(X5);
        X45 = diag(X4) - diag(X5);

        InteractionMtx = eta ./ sqrt(X12.^2 + Beta^2) + eta ./ sqrt(X13.^2 + Beta^2) + eta ./ sqrt(X14.^2 + Beta^2) + eta ./ sqrt(X15.^2 + Beta^2) + eta ./ sqrt(X23.^2 + Beta^2) + eta ./ sqrt(X24.^2 + Beta^2) + eta ./ sqrt(X25.^2 + Beta^2) + eta ./ sqrt(X34.^2 + Beta^2) + eta ./ sqrt(X35.^2 + Beta^2) + eta ./ sqrt(X45.^2 + Beta^2);
        InteractionMtx = sparse(1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5), 1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5), InteractionMtx);

        Hamiltonian = Hamiltonian + InteractionMtx;

        Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;

        Temp_Proj = zeros((Nx1 * Nx2 * Nx3 * Nx4 * Nx5), 3);


        k = 1;
        for j = 1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5)
            if (X12(j) < -eps) && (X23(j) < -eps) &&  (X34(j) < -eps) && (X45(j) < -eps)
                Temp_Proj(k, :) = [k, j, 1];
                k = k + 1;
            end
        end

        Temp_Proj = Temp_Proj(1:(k-1), :);
        Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, (Nx1 * Nx2 * Nx3 * Nx4 * Nx5));

        Hamiltonian = Projector * Hamiltonian * Projector';

        X1 = Projector * X1 * Projector';
        X2 = Projector * X2 * Projector';
        X3 = Projector * X3 * Projector';
        X4 = Projector * X4 * Projector';
        X5 = Projector * X5 * Projector';

        [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
        Spectrum    = diag(E);
        dE       = Spectrum(2) - Spectrum(1);

        Polarization.Splittings(AlphaInd, KappaInd) = dE;

        Int1 = Psi(:, 1)' * (X1 * Psi(:, 1));
        Int2 = Psi(:, 1)' * (X2 * Psi(:, 1));
        Int3 = Psi(:, 1)' * (X3 * Psi(:, 1));
        Int4 = Psi(:, 1)' * (X4 * Psi(:, 1));
        Int5 = Psi(:, 1)' * (X5 * Psi(:, 1));

        P = Int1 + Int2 + Int3 + Int4 + Int5;
        Polarization.Polarization(AlphaInd, KappaInd) = P;
        WaveFunction    = Projector' * Psi(:, 1);
        WFT = reshape(WaveFunction, [Nx1 Nx2 Nx3 Nx4 Nx5]);

        %         for Ind1 = 1:Nx1
        %             Psi1(Ind1) = 0;
        %             Psi5(Ind1) = 0;
        %             for Ind2 = 1: Nx2
        %                 Psi2(Ind2) = 0;
        %                 Psi4(Ind2) = 0;
        %                 for Ind3 = 1:Nx3
        %                     Psi3(Ind3) = 0;
        %                     for Ind4 = 1:Nx4
        %                         for Ind5 = 1:Nx5
        %                             Psi1(Ind1) = Psi1(Ind1) + WFT(Ind1, Ind2, Ind3, Ind4, Ind5);
        %                             Psi2(Ind2) = Psi2(Ind2) + WFT(Ind2, Ind1, Ind3, Ind4, Ind5);
        %                             Psi3(Ind3) = Psi3(Ind3) + WFT(Ind2, Ind3, Ind1, Ind4, Ind5);
        %                             Psi4(Ind4) = Psi4(Ind4) + WFT(Ind2, Ind3, Ind4, Ind1, Ind5);
        %                             Psi5(Ind5) = Psi5(Ind5) + WFT(Ind2, Ind3, Ind4, Ind5, Ind1);
        %                         end
        %                     end
        %                 end
        %             end
        %         end

        for Ind1 = 1:Nx1
            Psi1(Ind1) = 0;
            for Ind2 = 1: Nx2
                for Ind3 = 1:Nx3
                    for Ind4 = 1:Nx4
                        for Ind5 = 1:Nx5
                            Psi1(Ind1) = Psi1(Ind1) + WFT(Ind1, Ind2, Ind3, Ind4, Ind5);
                        end
                    end
                end
            end
        end

        for Ind2 = 1:Nx2
            Psi2(Ind2) = 0;
            for Ind1 = 1: Nx1
                for Ind3 = 1:Nx3
                    for Ind4 = 1:Nx4
                        for Ind5 = 1:Nx5
                            Psi2(Ind2) = Psi2(Ind2) + WFT(Ind1, Ind2, Ind3, Ind4, Ind5);
                        end
                    end
                end
            end
        end

        for Ind3 = 1:Nx3
            Psi3(Ind3) = 0;
            for Ind1 = 1: Nx1
                for Ind2 = 1:Nx2
                    for Ind4 = 1:Nx4
                        for Ind5 = 1:Nx5
                            Psi3(Ind3) = Psi3(Ind3) + WFT(Ind1, Ind2, Ind3, Ind4, Ind5);
                        end
                    end
                end
            end
        end

        for Ind4 = 1:Nx4
            Psi4(Ind4) = 0;
            for Ind1 = 1: Nx1
                for Ind2 = 1:Nx2
                    for Ind3 = 1:Nx3
                        for Ind5 = 1:Nx5
                            Psi4(Ind4) = Psi4(Ind4) + WFT(Ind1, Ind2, Ind3, Ind4, Ind5);
                        end
                    end
                end
            end
        end

        for Ind5 = 1:Nx5
            Psi5(Ind5) = 0;
            for Ind1 = 1: Nx1
                for Ind2 = 1:Nx2
                    for Ind3 = 1:Nx3
                        for Ind4 = 1:Nx4
                            Psi5(Ind5) = Psi5(Ind5) + WFT(Ind1, Ind2, Ind3, Ind4, Ind5);
                        end
                    end
                end
            end
        end

        Psi1 = Psi1 ./ norm(Psi1);
        Psi2 = Psi2 ./ norm(Psi2);
        Psi3 = Psi3 ./ norm(Psi3);
        Psi4 = Psi4 ./ norm(Psi4);
        Psi5 = Psi5 ./ norm(Psi5);

        figure(10)
        clf(figure(10))
        hold on
        mult = 30
        xxx = linspace(-10, 10, 100);
        plot(xxx, 10^-3 * (0.25 * (xxx.^2 - Alpha(AlphaInd)).^2 + Kappa(KappaInd) * xxx))
        plot(x1, abs(Psi1))
        plot(x2, abs(Psi2))
        plot(x3, abs(Psi3))
        plot(x4, abs(Psi4))
        plot(x5, abs(Psi5))

        xlim([-7, 7])
        ylim([-0.1 1])
        hold off
        Alpha(AlphaInd)
        Kappa(KappaInd)

        WF.WF1 = Psi1;
        WF.WF2 = Psi2;
        WF.WF3 = Psi3;
        WF.WF4 = Psi4;
        WF.WF5 = Psi5;
        WF.x1  = x1;
        WF.x2  = x2;
        WF.x3  = x3;
        WF.x4  = x4;
        WF.x5  = x5;
        WF.alpha    = Alpha(AlphaInd);
        WF.kappa    = Kappa(KappaInd);
        name = ['WaveFunction' num2str(10 * abs(Kappa(KappaInd)))];
        save(name, "WF")

    end
end


%%

figure(10)
clf(figure(10))
hold on
xxx = linspace(-10, 10, 100);
plot(xxx, 10^-3 * (0.25 * (xxx.^2 - Alpha(AlphaInd)).^2 + Kappa(KappaInd) * xxx))
plot(x1, abs(Psi1))
plot(x2, abs(Psi2))
plot(x3, abs(Psi3))
plot(x4, abs(Psi4))
plot(x5, abs(Psi5))
xlim([-7, 7])
ylim([-0.1 1])
hold off