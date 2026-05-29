function [TEMP, AL10N, AMMW, RHO] = cira(height_km, ALSA, DESA, ALSO, DESO, RJDAYS, GEO)
    % Константы и массивы
    AL10 = 2.3025851;
    ALPHA = [0.0, 0.0, 0.0, 0.0, -0.38];
    AMW = [28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797];
    AVOGAD = 6.02257e26;
    CONS25 = 0.35355339;
    FOURPI = 12.566371;
    TWOPI = 6.2831853;
    PIOV2 = 1.5707963;
    PIOV4 = 0.78539816;
    FRAC = [0.7811, 0.20955, 9.3432e-3, 6.1471e-6];
    RSTAR = 8314.32;
    R1 = 1.0e-2;
    R2 = 2.5e-2;
    R3 = 7.5e-2;
    WT = [0.31111111, 1.4222222, 0.53333333, 1.4222222, 0.31111111];

    % Перевод Юлианской даты
    AMJD = RJDAYS - 2400000.5;

    % Солнечные параметры
    TSUBC = 379.0 + 3.24 * GEO(2) + 1.3 * (GEO(1) - GEO(2));

    ETA = 0.5 * abs(DESA - DESO);
    THETA = 0.5 * abs(DESA + DESO);
    H = ALSA - ALSO;

    TAU = H - 0.64577182 + 0.10471976 * sin(H + 0.75049158);

    C = cos(ETA)^2.2;
    S = sin(THETA)^2.2;
    DF = S + (C - S) * abs(cos(0.5 * TAU))^3;
    TSUBL = TSUBC * (1.0 + 0.3 * DF);

    EXPKP = exp(GEO(3));
    DTG18 = 28.0 * GEO(3) + 0.03 * EXPKP;
    DTG20 = 14.0 * GEO(3) + 0.02 * EXPKP;
    DLR20 = 0.012 * GEO(3) + 1.2e-5 * EXPKP;

    F = 0.5 * (tanh(0.04 * (height_km - 350.0)) + 1.0);
    DLRGM = DLR20 * (1.0 - F);
    DTG = DTG20 * (1.0 - F) + DTG18 * F;
    TINF = TSUBL + DTG;

    TEMP = zeros(2, 1);
    TEMP(1) = TINF;  % Температура экзосферы

    % Промежуточные расчеты для TLOCAL
    TSUBX = 371.6678 + 5.18806e-2 * TINF - 294.3503 * exp(-2.16222e-3 * TINF);
    GSUBX = 0.054285714 * (TSUBX - 183.0);

    TC = zeros(4, 1);
    TC(1) = TSUBX;
    TC(2) = GSUBX;
    TC(3) = (TINF - TSUBX) / PIOV2;
    TC(4) = GSUBX / TC(3);

    % Заглушки на выходы (пока пустые)
    AL10N = zeros(6, 1);
    AMMW = NaN;
    RHO = NaN;

    % Следующий шаг: расчёт плотности и температур (от 90 км)

TC(4) = GSUBX / TC(3);

% --- Расчёт плотности от 90 до 100 км ---
    Z1 = 90.0;
    Z2 = min(height_km, 100.0);
    AL = log(Z2 / Z1);
    N = floor(AL / R1) + 1;
    ZR = exp(AL / N);

    AMBAR1 = ambar(Z1);
    TLOC1 = tlocal(Z1, TC);
    ZEND = Z1;
    SUM2 = 0.0;
    AIN = AMBAR1 * grav(Z1) / TLOC1;

    for i = 1:N
        Z = ZEND;
        ZEND = ZR * Z;
        DZ = 0.25 * (ZEND - Z);
        SUM1 = WT(1) * AIN;
        for j = 2:5
            Z = Z + DZ;
            AMBAR2 = ambar(Z);
            TLOC2 = tlocal(Z, TC);
            GRAVL = grav(Z);
            AIN = AMBAR2 * GRAVL / TLOC2;
            SUM1 = SUM1 + WT(j) * AIN;
        end
        SUM2 = SUM2 + DZ * SUM1;
    end

    FACT1 = 1e3 / RSTAR;
    RHO = 3.46e-6 * AMBAR2 * TLOC1 * exp(-FACT1 * SUM2) / (AMBAR1 * TLOC2);

    ANM = AVOGAD * RHO;
    AN = ANM / AMBAR2;

    FACT2 = ANM / 28.96;
    ALN = zeros(1, 6);
    ALN(1) = log(FRAC(1) * FACT2);
    ALN(4) = log(FRAC(3) * FACT2);
    ALN(5) = log(FRAC(4) * FACT2);
    ALN(2) = log(FACT2 * (1.0 + FRAC(2)) - AN);
    ALN(3) = log(2.0 * (AN - FACT2));

    if height_km <= 100.0
        TEMP(2) = TLOC2;
        ALN(6) = ALN(5) - 25.0;
        for i = 1:6
            AL10N(i) = ALN(i) / AL10;
        end
        AMMW = sum(exp(ALN) .* AMW) / sum(exp(ALN));
        RHO = sum(exp(ALN) .* AMW) / AVOGAD;
        return
    end
% --- Расчёт плотности от 100 до 500 км ---
    Z3 = min(height_km, 500.0);
    AL = log(Z3 / Z);
    N = floor(AL / R2) + 1;
    ZR = exp(AL / N);
    SUM2 = 0.0;
    AIN = GRAVL / TLOC2;

    for i = 1:N
        Z = ZEND;
        ZEND = ZR * Z;
        DZ = 0.25 * (ZEND - Z);
        SUM1 = WT(1) * AIN;
        for j = 2:5
            Z = Z + DZ;
            TLOC3 = tlocal(Z, TC);
            GRAVL = grav(Z);
            AIN = GRAVL / TLOC3;
            SUM1 = SUM1 + WT(j) * AIN;
        end
        SUM2 = SUM2 + DZ * SUM1;
    end

    % --- Если выше 500 км, дополнительный интеграл ---
    Z4 = max(height_km, 500.0);
    AL = log(Z4 / Z);
    R = R2;
    if height_km > 500.0
        R = R3;
    end
    N = floor(AL / R) + 1;
    ZR = exp(AL / N);
    SUM3 = 0.0;

    for i = 1:N
        Z = ZEND;
        ZEND = ZR * Z;
        DZ = 0.25 * (ZEND - Z);
        SUM1 = WT(1) * AIN;
        for j = 2:5
            Z = Z + DZ;
            TLOC4 = tlocal(Z, TC);
            GRAVL = grav(Z);
            AIN = GRAVL / TLOC4;
            SUM1 = SUM1 + WT(j) * AIN;
        end
        SUM3 = SUM3 + DZ * SUM1;
    end

    if height_km > 500.0
        T500 = TLOC3;
        TEMP(2) = TLOC4;
        ALTR = log(TLOC4 / TLOC2);
        FACT2 = FACT1 * (SUM2 + SUM3);
        HSIGN = -1.0;
    else
        T500 = TLOC4;
        TEMP(2) = TLOC3;
        ALTR = log(TLOC3 / TLOC2);
        FACT2 = FACT1 * SUM2;
        HSIGN = 1.0;
    end

    for i = 1:5
        ALN(i) = ALN(i) - (1 + ALPHA(i)) * ALTR - FACT2 * AMW(i);
    end

    AL10T5 = log10(T500);
    ALNH5 = (5.5 * AL10T5 - 39.4) * AL10T5 + 73.13;
    ALN(6) = AL10 * (ALNH5 + 6) + HSIGN * (log(TLOC4 / TLOC3) + FACT1 * SUM3 * AMW(6));
    CAPPHI = mod((AMJD - 36204.0) / 365.2422, 1.0);

    TAU = CAPPHI + 0.09544 * ((0.5 + 0.5 * sin(TWOPI * CAPPHI + 6.035))^1.65 - 0.5);
    GOFT = 0.02835 + 0.3817 * (1 + 0.4671 * sin(TWOPI * TAU + 4.137)) ...
                 * sin(FOURPI * TAU + 4.259);

    FOFZ = (5.876e-7 * height_km^2.331 + 6.328e-2) * exp(-2.868e-3 * height_km);
    DLRSA = FOFZ * GOFT;

    % Учет DLRSL (солнечно-ориентированная вариация)
    if height_km > 460.0
        DLRSL = 0.0;
    else
        DLRSL = 0.014 * (height_km - 90.0) * exp(-1.3e-3 * (height_km - 90)^2) ...
              * sign(DESA) * (sin(TWOPI * CAPPHI + 1.72) * sin(DESA))^2;
    end

    DLR = AL10 * (DLRGM + DLRSA + DLRSL);
    for i = 1:6
        ALN(i) = ALN(i) + DLR;
    end

    % Коррекция по гелию
    DLNHE = 0.65 * abs(DESO / 0.4091609) * ...
        ((sin(PIOV4 - 0.5 * DESA * sign(DESO)))^3 - CONS25);
    ALN(5) = ALN(5) + AL10 * DLNHE;

    % Вывод итоговых значений
    SUMN = 0.0;
    SUMNM = 0.0;
    for i = 1:6
        AN = exp(ALN(i));
        SUMN = SUMN + AN;
        SUMNM = SUMNM + AN * AMW(i);
        AL10N(i) = ALN(i) / AL10;
    end

    AMMW = SUMNM / SUMN;
    RHO = SUMNM / AVOGAD;
    end