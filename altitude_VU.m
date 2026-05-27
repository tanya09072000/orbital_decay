function altitude_VU()

current_file = fileparts(mfilename('fullpath'));

cx_path = fullfile(current_file, '..', 'calculation_Cx', 'matlab');

addpath(cx_path);

% Моделирование орбитального снижения КА.
% На каждом витке: интегрирование до восходящего узла, обновление SMA и высоты.

    params  = setup_params();
    data    = load_solar_data(params.data_file);
    results = run_simulation(params, data);
    plot_results(results, params.duration_days);
    print_summary(results);
end

% ================================================================
function params = setup_params()
    params.input_data    = [6908.888999, 0.001009, 97.401, 317.904, 357.882, 2.118];
    params.m_dry         = 306.4;    % масса КА без топлива, кг
    params.m_fuel        = 3.6;      % начальная масса топлива, кг
    params.m_total       = 306.4 + 3.6;  % полная стартовая масса, кг
   % params.Cx            = 3.56;     % коэффициент аэродинамического сопротивления
    params.Sm            = 2.25;     % площадь миделевого сечения, м²
    params.duration_days = 28.9;     % продолжительность расчёта, сут
    params.start_dt      = datetime(2025, 8, 31, 23, 18, 47);
    params.data_file     = 'Koeff_future.csv';
end

% ================================================================
function data = load_solar_data(filename)
    opts = detectImportOptions(filename, 'Encoding', 'UTF-8');
    opts.VariableNamingRule = 'preserve';
    opts = setvaropts(opts, 'DATE', 'InputFormat', 'yyyy-MM-dd');
    data = readtable(filename, opts);
end

% ================================================================
function results = run_simulation(params, data)
    C  = orb_constants();
    mu = C.mu;

    y_vec = vector_sostoyania(params.input_data);
    r0    = norm(y_vec(1:3));

    duration_s = params.duration_days * 86400;
    n_est      = ceil(duration_s / 5700) + 20;   % оценка числа витков + запас

    sma_arr    = zeros(1, n_est);
    height_arr = zeros(1, n_est);
    time_arr   = NaT(1, n_est);

    n            = 1;
    current_time = 0;
    current_dt   = params.start_dt;
    sma_arr(n)   = params.input_data(1);
    height_arr(n) = r0 - C.Re;
    time_arr(n)  = current_dt;

    % Предвычисление атмосферы

    h_grid = 100:1:1000;   % высоты, км

    rho_grid = zeros(size(h_grid));
    Cx_grid = zeros(size(h_grid));

    for i = 1:length(h_grid)

        [~, ~, rho_grid(i)] = cira_main( ...
        h_grid(i), ...
        dateshift(params.start_dt, 'start', 'day'), ...
        data);

        [Cx_grid(i), ~] = compute_cx( ...
        h_grid(i), ...
        0, ...
        dateshift(params.start_dt, 'start', 'day'), ...
        data);

    end


    while current_time < duration_s

        r = norm(y_vec(1:3));
        v = norm(y_vec(4:6));

        denom = 2/r - v^2/mu;
        if abs(denom) < 1e-12
            warning('run_simulation: орбита вырожденная (denom≈0), остановка на витке %d', n);
            break;
        end
        r_now   = 1 / denom;
        T_orbit = 2*pi * sqrt(r_now^3 / mu);

        % Плотность для текущей высоты и текущей даты
        %rho = cira_main(r - C.Re, dateshift(current_dt, 'start', 'day'), data);

        orbit_opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
                            'MaxStep', T_orbit/1000, 'Events', @ascendingNodeEvent);
        [~, ~, te, ye, ~] = ode89(...
            @(t, y) two_body(t, y, params.Sm, params.m_total, h_grid, rho_grid, Cx_grid), ...
            [0, T_orbit], y_vec, orbit_opts);

        if isempty(te)
            warning('run_simulation: восходящий узел не найден на витке %d', n);
            break;
        end

        y_vec        = ye(1, :)';
        current_time = current_time + te(1);
        current_dt   = params.start_dt + seconds(current_time);

        r = norm(y_vec(1:3));
        v = norm(y_vec(4:6));

        n = n + 1;
        sma_arr(n)    = 1 / (2/r - v^2/mu);
        height_arr(n) = r - C.Re;
        time_arr(n)   = current_dt;
    end

    results.sma    = sma_arr(1:n);
    results.height = height_arr(1:n);
    results.time   = time_arr(1:n);
end

% ================================================================
function plot_results(results, duration_days)
    figure;
    plot(results.time, results.sma, 'b-', 'LineWidth', 1.5);
    xlabel('Дата');
    xtickformat('dd-MMM-yyyy');
    ylabel('Большая полуось, км');
    title(sprintf('Изменение большой полуоси за %.1f сут', duration_days));
    grid on;
end

% ================================================================
function print_summary(results)
    dh = results.sma(1) - results.sma(end);
    fprintf('Изменение SMA:   %.6f км\n', dh);
    fprintf('Период расчёта:  %s — %s\n', ...
        string(results.time(1),   'yyyy-MM-dd HH:mm:ss'), ...
        string(results.time(end), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('Начальная SMA:   %.6f км\n', results.sma(1));
    fprintf('Конечная SMA:    %.6f км\n', results.sma(end));

    T = table(results.time(:), results.sma(:), results.height(:), ...
        'VariableNames', {'Дата', 'Большая_полуось_км', 'Высота_км'});
end
