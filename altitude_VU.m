function altitude_VU()
% программа показывает насколько уменьшилась большая полуось в ВУ
% использует доп. функции по поиску ВУ
% заменить шаг в один орбитальный период на поиск восходящего узла как события
input_data = [6908.888999 0.001009 97.401 317.904 357.882 2.118];
Re = 6378.137; % радиус Земли на экваторе, км
mu = 398600.4415; % гравитационный параметр Земли, км³/с²
g0 = 9.80665; % ускорение свободного падения, м/с²
%Cx = 4.25; % коэффициент аэродинамического сопротивления
%Sm = 1.87; % площадь сечения, м²
m_KA = 310; % масса КА, кг
Cx = 3.56; % коэф. аэрод. сопротивления
Sm = 2.25;  % площадь сечения, м² 
mass_fuel = 3.6; % масса топлива, кг
m_KA = m_KA - mass_fuel;

duration_days = 28.9; % кол-во дней 6.9/
duration_seconds = duration_days * 86400; 

y_vec = vector_sostoyania(input_data); % начальный вектор-состояния КА [x; y; z; vx; vy; vz]
HEGEO = norm(y_vec(1:3)) - Re; % высота орбиты, км

% Инициализация
m_total = m_KA + mass_fuel;
start_datetime = datetime(2025, 08, 31, 23, 18, 47); % Инициализация текущей даты (год, месяц, день, часы, минуты, секунды)
current_datetime = start_datetime;
current_time = 0;

semi_major_axis = [];
height = []; % высота, км
modRV = []; % модуль вектора состояния
time = datetime.empty; % время, сек
Y_vec = []; % массив для векторов состояния
crossing_times = datetime.empty;

semi_major_axis(end+1) = input_data(1);
height(end+1) = HEGEO;
time(end+1) = current_datetime;
modRV(end+1) = norm(y_vec);
crossing_times(end+1) = current_datetime;

% Формирование таблицы плотностей по месяцам
startDate = datetime('2025-09-01', 'InputFormat', 'yyyy-MM-dd');
endDate = datetime('2040-01-01', 'InputFormat', 'yyyy-MM-dd');
opts = detectImportOptions('Koeff_future.csv', 'Encoding', 'UTF-8');
opts.VariableNamingRule = 'preserve';
opts = setvaropts(opts, 'DATE', 'InputFormat', 'yyyy-MM-dd');
data = readtable('Koeff_future.csv', opts);
allDate = startDate:calmonths(1):endDate;
allRHO = zeros(length(allDate), 1);
for k = 1:length(allDate)
    [~, ~, allRHO(k)] = cira_main(HEGEO, allDate(k), data);
end

while current_time < duration_seconds

    % Определить плотность для текущей даты
    [~, date_idx] = min(abs(allDate - dateshift(current_datetime, 'start', 'month')));
    rho = allRHO(date_idx);

    % Расчет орбитального периода
    denom = 2/norm(y_vec(1:3)) - norm(y_vec(4:6))^2 / mu;
    if abs(denom) < 1e-12
        warning('Знаменатель vis-viva близок к нулю — орбита параболическая или гиперболическая. Прерывание.');
        break;
    end
    r_now = 1 / denom;
    T_orbit = 2 * pi * sqrt(r_now^3 / mu);
  
    % Интегрирование движения на один виток
    tspan = [0, T_orbit];
    options = odeset('RelTol',1e-10, 'AbsTol',1e-12, 'MaxStep', T_orbit/1000, 'Events', @ascendingNodeEvent);
    [t, y, te, ye, ie] = ode89(@(t, y) two_body(t, y, Cx, Sm, m_total, rho), tspan, y_vec, options);

    if isempty(te)
        warning('Восходящий узел не найден, прерывание интегрирования');
        break;
    end

    % Обновление состояния
    t_node = te(1);          % время наступления события
    y_vec = ye(1, :)';       % состояние в момент узла
    crossing_times(end+1) = start_datetime + seconds(current_time + t_node);

    current_time = current_time + t_node; % время в точке восходящего узла
    current_datetime = start_datetime + seconds(current_time);

    % Сохраняем высоту и время
    denom_a = 2/norm(y_vec(1:3)) - norm(y_vec(4:6))^2 / mu;
    if abs(denom_a) < 1e-12
        warning('Знаменатель vis-viva близок к нулю при сохранении SMA. Прерывание.');
        break;
    end
    a = 1 / denom_a; % большая полуось, км
    semi_major_axis(end+1) = a;
    height(end+1) = norm(y_vec(1:3)) - Re; % высота в восходящем узле
    time(end+1) = current_datetime;

    modRV(end+1) = norm(y_vec); 
end

dh = semi_major_axis(1)-semi_major_axis(end); 

% График
figure;
plot(time, semi_major_axis, 'b-', 'LineWidth', 1.5);
xlabel('Дата');
datetick('x','dd-mmm-yyyy','keepticks'); % формат даты на оси X
ylabel('SMA');
title(sprintf('Изменение SMA орбиты за %d дней', duration_days));
grid on;

fprintf('Высота упала на %.6f км\n', dh);
fprintf('Период расчёта: с %s по %s\n',datestr(time(1), 'yyyy-mm-dd HH:MM:SS'), datestr(time(end), 'yyyy-mm-dd HH:MM:SS'));
fprintf('Начальная SMA %.6f км\n', semi_major_axis(1));
fprintf ('Конечная SMA %.6f км\n', semi_major_axis(end));
T = table(crossing_times(:), semi_major_axis(:), height(:), 'VariableNames', {'Дата', 'Большая полуось_км', 'Высота_км'});
% disp(T)

end