function [TEMP, AMMW, RHO,f107,f81,Kp] = cira_main(height_km, date_str)

        % Загрузка данных
        koeffPath = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'Koeff.csv');
        opts = detectImportOptions(koeffPath, 'Encoding', 'UTF-8');
        opts.VariableNamingRule = 'preserve';
        opts = setvaropts(opts, 'DATE', 'InputFormat', 'yyyy-MM-dd');
        data = readtable(koeffPath, opts);

        try
            myDate = datetime(date_str, 'InputFormat', 'yyyy-MM-dd');
        catch
            error('Неверный формат даты. Введите, например: 2020-01-04');
        end

        % --- Поиск F10.7 по дате ---
        idx = data.DATE == myDate;
        if ~any(idx)
            error('Дата не найдена в таблице солнечных данных.');
        end

        f107 = data.("F10.7_OBS")(idx);  % колонка с F10.7 по наблюдениям
        f81 = data.("F10.7_OBS_LAST81")(idx); % колонка с F10.7 по наблюдениям 81 дня
        Kp = data.("KP_SUM")(idx);
        Kp = Kp/80;
        GEO = [f107; f81; Kp];  % GEO(3) — геомагнитный индекс, по умолчанию 1.0
        %GEO = [150; 150; 3];
        %disp (GEO);
        ALSA = 0.0;
        DESA = 0.0;
        ALSO = 0.0;
        DESO = pi/2;
        RJDAYS = juliandate(myDate);
        

        TEMP = zeros(2, 1);
        AL10N = zeros(6, 1);

        [TEMP, AL10N, AMMW, RHO] = cira(height_km, ALSA, DESA, ALSO, DESO, RJDAYS, GEO);

        AMMW = AMMW;
        %fprintf('Температура экзосферы и локальная температура T, K: %.2f %.2f\n', TEMP(1), TEMP(2));
        %fprintf('Молекулярная масса: %.6f\n', AMMW);
        %fprintf('Плотность, кг/м^3: %.6e\n', RHO);
end