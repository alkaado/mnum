function main()
    %n = 1000;
    %dt = 0.006;
    %vx = linspace(0, dt, n);
    %vy = test2(vx);
    %X = DFT(vy, n);
    n = 1024;
    timeframe = 1;
    vx = linspace(0, timeframe, n);
    vy = test2(vx);
    X = abs(FFT(vy, n));
    plot(fftfreq(n, timeframe/n), 2 * abs(X/n));
    %X2 = fft(vy,n);
    %plot(linspace(0, round(n/8)/dt, round(n/8)), 2 .* abs(X(1:round(n/8))./n))
end

function vf = fftfreq(n, dt)
    tmp = [];
    
    for i = 1: n / 2
           tmp(i) = i / (n * dt);
           tmp(i + n/2) = (- n/2 + i) / (n * dt);
    end
    
    vf = tmp;
end

function vy = test2(vx)
    vy = 512 * sin(2 * pi * 80 .* vx);
end

function vy = test(vx)
    vy = 5 * sin((2 * pi * 200 ) .* vx) + (10 * sin((2 * pi * 20000 ) .* vx));
end

function X = DFT(xt, n)
    a0 = 0;
    a1 = 0;
    a2 = 0;
    X = zeros(n,1);
    S = 0;
    a0 = 2 * pi * 1/n;
    for k = 1:n
        a1 = a0 * (k - 1);
        for i = 1:n
            a2 = a1 * (i - 1);
            
            S = S + xt(i) * cos(a2) + -xt(i) * sin(a2) * 1j;
            %fprintf("%d %f %f\n", i, real(S), imag(S));
        end

        fprintf("%d %f %f %f \n", k, real(S), imag(S), abs(S));
        X(k) = S;
        S = 0;
    end
end

function X = FFT(x, n)
    WN = exp(-j * 2 * pi / n);
    N = round(n / 2);
    S1 = 0j;
    S2 = 0j;
    WNK = WN;
    WNK2 = 1j;
    WNK2I = 1j;
    
    for k = 1:n
        for i = 1:N
            S1 = S1 + x((i*2) - 1) * WNK2I;
            S2 = S2 + x(i*2) * WNK2I;
            WNK2I = WNK2I * WNK2;
        end
        X(k) = S1 + S2 * WNK;
        S1 = 0;
        S2 = 0;
        WNK = WNK * WN;
        WNK2 = WNK * WNK;
    end
end

function X = FFT2(x, n)
    WN = exp(-j * 2 * pi / n);
    N = round(n / 2);
    S = 0j;
    WNK = 1j;
    WNK2 = 1j;
    
    for k = 1:n
        WNK = WN^k;
        WNK2 = WNK^2;
        WNK2I = 1;
        for i = 1:N
            S = S + (x((i*2) - 1) + x(i*2) * WNK) * WNK2I; 
            WNK2I = WNK2I * WNK2;
        end
        X(k) = S;
        S = 0;
        
    end
end
