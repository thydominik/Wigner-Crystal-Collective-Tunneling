clc
clear all
ED = 0;
disp('WavePacket simulation')

NoP = 2^11;
alpha = 20;

x               = linspace(-10, 10, NoP);
IntervalLength  = abs(max(x) - min(x));

dx = x(2) - x(1);

k       = zeros(1, NoP);
dk      = 2 * pi/IntervalLength;
Omega   = dk * NoP;
k       = linspace(-Omega/2, Omega/2, NoP);

dt = 10^-3;

k0      = 0.0;
sigma   = 2;
x0      = -sqrt(alpha);

V = 0.25 * (x.^2 - alpha).^2;

ImpulseOp = fftshift(exp(-1i*dt*k.^2 / 4));

IterSteps = 2*10^3;
for j = 1:IterSteps
    if j == 1
        if ED == 1
            PotentialMtx = sparse(diag(V));
            for i = 2:NoP
                K = 1/(2 * dx^2);
                KineticMtx(i - 1, i) = -K;
                KineticMtx(i, i - 1) = -K;
            end
                Hamiltonian = PotentialMtx + sparse(KineticMtx);
    
                [Psi, E]    = eig(full(Hamiltonian));
                E           = diag(E);
                psi0 = Psi(:, 1);
                figure(1)
                clf(figure(1))
                hold on
                    plot(x, psi0)
                hold off
        else
            psi0= exp((x - x0).^2 / -(sigma^2)) .* exp(1i * k0 * x);
            psi0 = psi0 ./ norm(psi0);
        end
    end

    psi = fft(psi0);

    for h = 1:NoP
        psi(h) = ImpulseOp(h) * psi(h);
    end
    psi = ifft(psi);
    for h = 1:NoP
        psi(h) = exp(-1i * dt * V(h)) * psi(h);
    end
    psi = fft(psi);
    for h = 1:NoP
        psi(h) = ImpulseOp(h) * psi(h);
    end

    psi = ifft(psi);
    psi0 = psi;
    figure(2)
    clf(figure(2))
    hold on
        plot(x, V ./ norm(V))
        plot(x, abs(psi ./ norm(psi)))
        
    hold off
    T(j) = sum(abs(psi(end/2 : end)).^2);
end
%%
figure(3)
clf(figure(3))
hold on
%set(gca, 'Yscale', 'log')
plot(linspace(0, IterSteps * dt, IterSteps), T)
hold off