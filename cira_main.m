function [TEMP, AMMW, RHO,f107,f81,Kp] = cira_main(HEGEO, userDate, data)

        % Поиск F10.7 по дате (data уже загружена)
        d = dateshift(userDate, 'start', 'day');
        idx = year(data.DATE)  == year(d)  & ...
        month(data.DATE) == month(d) & ...
        day(data.DATE)   == day(d);
            if ~any(idx)
                error('Дата не найдена в таблице солнечных данных: %s', string(userDate));
            end
    
        f107 = data.("F10.7_OBS")(idx);  % колонка с F10.7 по наблюдениям
        f81 = data.("F10.7_OBS_LAST81")(idx); % колонка с F10.7 по наблюдениям 81 дня
        Kp = data.("KP_SUM")(idx);
        Kp = Kp/80;
        GEO = [f107; f81; Kp];  % GEO(3) — геомагнитный индекс, по умолчанию 1.0
        %GEO = [150; 150; 3];
        %disp (GEO);
        
        day_of_year = day(userDate, 'dayofyear');
        DESA = deg2rad(23.44 * sin(2*pi * (day_of_year - 80) / 365));  % солнечное склонение
        
        ALSA = 0.0;
        ALSO = 0.0;
        DESO = pi/2;
        RJDAYS = juliandate(userDate);
        

        TEMP = zeros(2, 1);
        AL10N = zeros(6, 1);

        [TEMP, AL10N, AMMW, RHO] = cira(HEGEO, ALSA, DESA, ALSO, DESO, RJDAYS, GEO);

        AMMW = AMMW;
        %fprintf('Температура экзосферы и локальная температура T, K: %.2f %.2f\n', TEMP(1), TEMP(2));
        %fprintf('Молекулярная масса: %.6f\n', AMMW);
        %fprintf('Плотность, кг/м^3: %.6e\n', RHO);
end