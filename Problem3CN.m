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
    dt = dx;
    dt = Tmax / (ceil(Tmax / dt));
    
    xvals = ((0:N) * dx) + xmin;
    tvals = (0:(Tmax / dt)) * dt;
    
    [T, X] = meshgrid(tvals, xvals);
    
    U = zeros(length(xvals)-2, length(tvals));
    U(:, 1) = u0(xvals(2:end-1));
    
    mainDiag = 1 + (dt/(dx^2)) * ones(1, N-1);
    subDiags = -0.5 * (dt / (dx^2)) * ones(1, N-2);
    A = diag(mainDiag, 0) + diag(subDiags, -1) + diag(subDiags, 1);
    
    for j = 2:length(tvals)
        Uprev = (1 - (dt / (dx^2))) * U(:, j-1) + 0.5 * (dt / (dx^2)) * [U(2:end, j-1); 0] + 0.5 * (dt / (dx^2)) * [0; U(1:end-1, j-1)];
        Uprev(1) = Uprev(1) + 0.5 * (dt / (dx^2)) * uNegBdry(tvals(j)) + 0.5 * (dt / (dx^2)) * uNegBdry(tvals(j-1));
        Uprev(end) = Uprev(end) + 0.5 * (dt / (dx^2)) * uPosBdry(tvals(j)) + 0.5 * (dt / (dx^2)) * uPosBdry(tvals(j-1));
        
        U(:, j) = A \ Uprev;
    end
    
    Ufull = [uNegBdry(tvals); U; uPosBdry(tvals)];
    
    trueSol = v(T+1, X);
    error = abs(trueSol(:, end) - Ufull(:, end));
    L2errors(iter) = sqrt(dx)*norm(error, 2);
    LinfErrors(iter) = norm(error, "inf");

end

dxvals = (xmax - xmin) ./ Nvals;

L2coeff = polyfit(log(dxvals), log(L2errors), 1);
LinfCoeff = polyfit(log(dxvals), log(LinfErrors), 1);
ErrorOrderForL2NormCN = L2coeff(1)
ErrorOrderForLinfNormCN = LinfCoeff(1)