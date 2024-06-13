xmin = -5;
xmax = 5;
Tmax = 3;
v = @(t, x) (1 ./ sqrt(4 * pi * t)) .* exp(-((x-2).^2) ./ (4 * t));

u0 = @(x) v(1, x);
uNegBdry = @(t) v(t + 1, -5);
uPosBdry = @(t) v(t + 1, 5);

Nvals = [20 40 80 160];
L2errors = zeros(1, length(Nvals));
LinfErrors = zeros(1, length(Nvals));

for iter = 1:length(Nvals)
    N = Nvals(iter);
    dx = (xmax - xmin) / N;
    dt = 0.4 * dx * dx;
    dt = Tmax / (ceil(Tmax / dt));
    
    xvals = ((0:N) * dx) + xmin;
    tvals = (0:(Tmax / dt)) * dt;
    [T, X] = meshgrid(tvals, xvals);
    
    U = zeros(length(xvals), length(tvals));
    U(:, 1) = u0(xvals);
    U(1, :) = uNegBdry(tvals);
    U(end, :) = uPosBdry(tvals);
    
    for k = 2:length(tvals)
        for j = 2:(length(xvals)-1)
            U(j, k) = (dt / (dx^2)) * (U(j+1, k-1) - 2*U(j, k-1) + U(j-1, k-1)) + U(j, k-1);
        end
    end
    
    trueSol = v(T+1, X);
    error = trueSol(:, end) - U(:, end);
    L2errors(iter) = sqrt(dx)*norm(error, 2);
    LinfErrors(iter) = norm(error, "inf");
end
dxvals = (xmax - xmin) ./ Nvals;
L2coeff = polyfit(log(dxvals), log(L2errors), 1);
LinfCoeff = polyfit(log(dxvals), log(LinfErrors), 1);
ErrorOrderForL2NormExplicitEuler = L2coeff(1)
ErrorOrderForLinfNormExplicitEuler = LinfCoeff(1)