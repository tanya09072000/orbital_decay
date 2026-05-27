function [TEMP, AL10N, AMMW, RHO] = cira(alt_km, ALSA, DESA, ALSO, DESO, jd, GEO)
% Модель атмосферы CIRA-86.
% Вход:  alt_km — геоцентрическая высота, км (>= 90)
%        ALSA   — прямое восхождение Солнца, рад
%        DESA   — склонение Солнца, рад
%        ALSO   — прямое восхождение орбитального узла, рад
%        DESO   — склонение орбитального узла, рад
%        jd     — юлианская дата
%        GEO    — [F10.7; F10.7_81day; Kp]
% Выход: TEMP   — [T_inf; T_local], К
%        AL10N  — log10 числовых концентраций шести компонент
%        AMMW   — средняя молярная масса, г/моль
%        RHO    — плотность, кг/м³

    % --- Модельные константы ---
    AL10   = log(10);
    ALPHA  = [0.0, 0.0, 0.0, 0.0, -0.38];
    AMW    = [28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797];
    AVOGAD = 6.02257e26;
    CONS25 = 0.35355339;
    FRAC   = [0.7811, 0.20955, 9.3432e-3, 6.1471e-6];
    RSTAR  = 8314.32;
    PIOV2  = pi / 2;
    PIOV4  = pi / 4;
    TWOPI  = 2 * pi;
    FOURPI = 4 * pi;
    R1 = 1.0e-2;   % шаг интегрирования log(z) для слоя  90-100 км
    R2 = 2.5e-2;   % шаг для слоя 100-500 км
    R3 = 7.5e-2;   % шаг для слоя > 500 км

    % --- Температурная модель ---
    AMJD  = jd - 2400000.5;
    TSUBC = 379.0 + 3.24*GEO(2) + 1.3*(GEO(1) - GEO(2));

    ETA   = 0.5 * abs(DESA - DESO);
    THETA = 0.5 * abs(DESA + DESO);
    H     = ALSA - ALSO;
    TAU   = H - 0.64577182 + 0.10471976 * sin(H + 0.75049158);
    C_ang = cos(ETA)^2.2;
    S_ang = sin(THETA)^2.2;
    DF    = S_ang + (C_ang - S_ang) * abs(cos(0.5*TAU))^3;
    TSUBL = TSUBC * (1.0 + 0.3*DF);

    EXPKP = exp(GEO(3));
    DTG18 = 28.0*GEO(3) + 0.03*EXPKP;
    DTG20 = 14.0*GEO(3) + 0.02*EXPKP;
    DLR20 = 0.012*GEO(3) + 1.2e-5*EXPKP;

    F     = 0.5 * (tanh(0.04*(alt_km - 350.0)) + 1.0);
    DLRGM = DLR20 * (1.0 - F);
    DTG   = DTG20*(1.0 - F) + DTG18*F;
    T_inf = TSUBL + DTG;

    TEMP    = zeros(2, 1);
    TEMP(1) = T_inf;

    T_sub = 371.6678 + 5.18806e-2*T_inf - 294.3503*exp(-2.16222e-3*T_inf);
    G_sub = 0.054285714 * (T_sub - 183.0);

    TC    = zeros(4, 1);
    TC(1) = T_sub;
    TC(2) = G_sub;
    TC(3) = (T_inf - T_sub) / PIOV2;
    if abs(TC(3)) < 1e-6
        TC(4) = 0.0;
    else
        TC(4) = G_sub / TC(3);
    end

    AL10N = zeros(6, 1);
    AMMW  = NaN;
    RHO   = NaN;

    % --- Интегрирование слоя 90-100 км ---
    Z1     = 90.0;
    Z2     = min(alt_km, 100.0);
    AMBAR1 = ambar(Z1);
    T_loc1 = tlocal(Z1, TC);
    AIN0   = AMBAR1 * grav(Z1) / T_loc1;
    [SUM_90_100, ~, G2, T2] = boole_integrate(Z1, Z2, AIN0, TC, R1, true);
    AMBAR2 = ambar(Z2);

    FACT1 = 1e3 / RSTAR;
    RHO   = 3.46e-6 * AMBAR2 * T_loc1 * exp(-FACT1*SUM_90_100) / (AMBAR1 * T2);

    ANM   = AVOGAD * RHO;
    AN    = ANM / AMBAR2;
    FACT2 = ANM / 28.96;

    ALN    = zeros(1, 6);
    ALN(1) = log(FRAC(1) * FACT2);
    ALN(4) = log(FRAC(3) * FACT2);
    ALN(5) = log(FRAC(4) * FACT2);
    ALN(2) = log(FACT2*(1.0 + FRAC(2)) - AN);
    ALN(3) = log(2.0*(AN - FACT2));

    if alt_km <= 100.0
        TEMP(2) = T2;
        ALN(6)  = ALN(5) - 25.0;
        AL10N   = (ALN / AL10)';
        exp_ALN = exp(ALN);
        AMMW    = sum(exp_ALN .* AMW) / sum(exp_ALN);
        RHO     = sum(exp_ALN .* AMW) / AVOGAD;
        return
    end

    % --- Интегрирование слоя 100-500 км ---
    Z3 = min(alt_km, 500.0);
    [SUM_100_500, ~, G3, T3] = boole_integrate(Z2, Z3, G2/T2, TC, R2, false);

    % --- Интегрирование слоя 500+ км (или от alt до 500 если alt < 500) ---
    Z4 = max(alt_km, 500.0);
    R  = R2;
    if alt_km > 500.0, R = R3; end
    [SUM_500p, ~, ~, T4] = boole_integrate(Z3, Z4, G3/T3, TC, R, false);

    if alt_km > 500.0
        T500    = T3;
        TEMP(2) = T4;
        ALTR    = log(T4 / T2);
        FACT2   = FACT1 * (SUM_100_500 + SUM_500p);
        HSIGN   = -1.0;
    else
        T500    = T4;
        TEMP(2) = T3;
        ALTR    = log(T3 / T2);
        FACT2   = FACT1 * SUM_100_500;
        HSIGN   = 1.0;
    end

    for k = 1:5
        ALN(k) = ALN(k) - (1 + ALPHA(k))*ALTR - FACT2*AMW(k);
    end

    AL10T5 = log10(T500);
    ALNH5  = (5.5*AL10T5 - 39.4)*AL10T5 + 73.13;
    ALN(6) = AL10*(ALNH5 + 6) + HSIGN*(log(T4/T3) + FACT1*SUM_500p*AMW(6));

    % --- Солнечные и сезонные поправки ---
    CAPPHI = mod((AMJD - 36204.0) / 365.2422, 1.0);
    TAU    = CAPPHI + 0.09544 * ((0.5 + 0.5*sin(TWOPI*CAPPHI + 6.035))^1.65 - 0.5);
    GOFT   = 0.02835 + 0.3817*(1 + 0.4671*sin(TWOPI*TAU + 4.137)) * sin(FOURPI*TAU + 4.259);
    FOFZ   = (5.876e-7*alt_km^2.331 + 6.328e-2) * exp(-2.868e-3*alt_km);
    DLRSA  = FOFZ * GOFT;

    if alt_km > 460.0
        DLRSL = 0.0;
    else
        % sign(DESA)*(sin(...)*sin(DESA))^2 переписано как sin(DESA)*|sin(DESA)|*sin(...)^2
        % — эквивалентно оригиналу, но непрерывно при DESA = 0
        DLRSL = 0.014*(alt_km - 90.0)*exp(-1.3e-3*(alt_km - 90)^2) ...
              * sin(TWOPI*CAPPHI + 1.72)^2 * sin(DESA)*abs(sin(DESA));
    end

    DLR  = AL10 * (DLRGM + DLRSA + DLRSL);
    ALN  = ALN + DLR;

    DLNHE  = 0.65*abs(DESO/0.4091609) * ((sin(PIOV4 - 0.5*DESA*sign(DESO)))^3 - CONS25);
    ALN(5) = ALN(5) + AL10*DLNHE;

    % --- Итоговая плотность ---
    exp_ALN = exp(ALN);
    SUMN    = sum(exp_ALN);
    SUMNM   = sum(exp_ALN .* AMW);
    AL10N   = (ALN / AL10)';

    if SUMN == 0
        warning('cira: нулевая суммарная концентрация на высоте %.0f км (underflow)', alt_km);
        AMMW = NaN;
    else
        AMMW = SUMNM / SUMN;
    end
    RHO = SUMNM / AVOGAD;
end

% -----------------------------------------------------------------------
function [integral, ZEND_out, G_out, T_out] = boole_integrate(Z_from, Z_to, AIN_init, TC, step_size, include_ambar)
% Интегрирование по логарифмической сетке высот методом Boole (5-точечные веса).
% Вход:  Z_from, Z_to    — диапазон высот, км
%        AIN_init        — значение подынтегрального выражения в Z_from
%        TC              — коэффициенты температурной модели
%        step_size       — максимальный шаг по log(z) (R1/R2/R3)
%        include_ambar   — true для слоя 90-100 км (включает молекулярную массу)
% Выход: integral        — значение интеграла
%        ZEND_out        — конечная высота (= Z_to)
%        G_out, T_out    — grav и tlocal в конечной точке

    WT       = [0.31111111, 1.4222222, 0.53333333, 1.4222222, 0.31111111];
    N        = floor(log(Z_to / Z_from) / step_size) + 1;
    ZR       = exp(log(Z_to / Z_from) / N);
    integral = 0.0;
    AIN      = AIN_init;
    ZEND_out = Z_from;
    G_out    = grav(Z_from);
    T_out    = tlocal(Z_from, TC);

    for k = 1:N
        Z_seg    = ZEND_out;
        ZEND_out = ZR * Z_seg;
        DZ       = 0.25 * (ZEND_out - Z_seg);
        S        = WT(1) * AIN;
        Z        = Z_seg;
        for j = 2:5
            Z     = Z + DZ;
            T_out = tlocal(Z, TC);
            G_out = grav(Z);
            if include_ambar
                AIN = ambar(Z) * G_out / T_out;
            else
                AIN = G_out / T_out;
            end
            S = S + WT(j) * AIN;
        end
        integral = integral + DZ * S;
    end
end
