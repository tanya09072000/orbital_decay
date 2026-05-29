function [TEMP, AMMW, RHO,f107,f81,Kp] = cira_main(HEGEO, userDate, data, lon_rad)
% lon_rad — географическая долгота КА, рад (необязательный аргумент)

        if nargin < 4 || isempty(lon_rad)
            lon_rad = 0;
        end

        % Поиск F10.7 по дате (data уже загружена)
        d = dateshift(userDate, 'start', 'day');
        idx = year(data.DATE)  == year(d)  & ...
        month(data.DATE) == month(d) & ...
        day(data.DATE)   == day(d);
            if ~any(idx)
                error('Дата не найдена в таблице солнечных данных: %s', string(userDate));
            end

        f107 = data.("F10.7_OBS")(idx);
        f81  = data.("F10.7_OBS_LAST81")(idx);
        Kp   = data.("KP_SUM")(idx) / 8;   % среднесуточный Kp (8 отсчётов/сут)
        GEO  = [f107; f81; Kp];

        day_of_year = day(userDate, 'dayofyear');
        DESA = deg2rad(23.44 * sin(2*pi * (day_of_year - 80) / 365));

        RJDAYS = juliandate(userDate);

        % --- Прямое восхождение Солнца (точность ~1°) ---
        T_jc    = (RJDAYS - 2451545.0) / 36525.0;
        M_sun   = deg2rad(mod(357.52911 + 35999.05029 * T_jc, 360));
        L_sun   = mod(280.46646 + 36000.76983 * T_jc, 360);
        lambda  = deg2rad(L_sun + 1.914602*sin(M_sun) + 0.019993*sin(2*M_sun));
        epsilon = deg2rad(23.439 - 0.013 * T_jc);
        ALSO    = atan2(cos(epsilon) * sin(lambda), cos(lambda));   % RA☉, рад

        % --- Гринвичское звёздное время → местное звёздное время КА ---
        GMST_rad = deg2rad(mod(280.46061837 + 360.98564736629 * (RJDAYS - 2451545.0), 360));
        ALSA     = mod(GMST_rad + lon_rad, 2*pi);

        DESO = pi/2;
        

        TEMP = zeros(2, 1);
        AL10N = zeros(6, 1);

        [TEMP, AL10N, AMMW, RHO] = cira(HEGEO, ALSA, DESA, ALSO, DESO, RJDAYS, GEO);

        AMMW = AMMW;
        %fprintf('Температура экзосферы и локальная температура T, K: %.2f %.2f\n', TEMP(1), TEMP(2));
        %fprintf('Молекулярная масса: %.6f\n', AMMW);
        %fprintf('Плотность, кг/м^3: %.6e\n', RHO);
end