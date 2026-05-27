function [TEMP, AMMW, RHO] = cira_main(alt_km, date, data)
% Обёртка над моделью CIRA-86: температура, молярная масса и плотность атмосферы.
% Вход:  alt_km — геоцентрическая высота, км (>= 90)
%        date   — datetime или строка 'yyyy-MM-dd' с датой расчёта
%        data   — таблица солнечных данных (из Koeff_future.csv)
% Выход: TEMP   — [T_inf; T_local], К
%        AMMW   — средняя молярная масса, г/моль
%        RHO    — плотность, кг/м³

    if isa(date, 'datetime')
        myDate = dateshift(date, 'start', 'day');
    elseif ischar(date) || isstring(date)
        myDate = datetime(date, 'InputFormat', 'yyyy-MM-dd');
    else
        error('cira_main: аргумент date должен быть datetime или строкой.');
    end

    idx = data.DATE == myDate;
    if ~any(idx)
        error('cira_main: дата %s не найдена в таблице солнечных данных.', ...
              string(myDate, 'yyyy-MM-dd'));
    end

    f107 = data.("F10.7_OBS")(idx);
    f81  = data.("F10.7_OBS_LAST81")(idx);
    % KP_SUM — сумма 8 трёхчасовых значений Kp (в единицах 0.1 Kp);
    % делим на 80 для перевода в средний Kp стандартной шкалы (0-9)
    Kp   = data.("KP_SUM")(idx) / 80;

    GEO = [f107; f81; Kp];

    % Приближённое склонение Солнца (точность ±1°) для учёта сезонной вариации
    doy  = day(myDate, 'dayofyear');
    DESA = deg2rad(23.45 * sin(2*pi * (284 + doy) / 365));

    ALSA  = 0.0;    % местное солнечное время не задано — усреднение по долготе
    ALSO  = 0.0;
    DESO  = pi / 2;
    jd    = juliandate(myDate);

    [TEMP, ~, AMMW, RHO] = cira(alt_km, ALSA, DESA, ALSO, DESO, jd, GEO);
end
