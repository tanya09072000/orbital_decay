function EO = element(Y1)
    % Преобразует вектор состояния в орбитальные элементы
    % Вход: input_data — [x, y, z, Vx, Vy, Vz]
    % Выход: EO — [a, e, i, om, w, nu]
    
    % Извлекаем компоненты
    x = Y1(1); y = Y1(2); z = Y1(3);
    Vx = Y1(4); Vy = Y1(5); Vz = Y1(6);
    
    r_vec = vpa([x; y; z]);
    v_vec = vpa([Vx; Vy; Vz]);
    mu = double(398600.4415);  % гравитационный параметр Земли (км^3/с^2)

    % Модуль вектора положения и скорости
    r = norm(r_vec);
    v = norm(v_vec);

    % Орбитальная энергия и большая полуось
    energy = v^2 / 2 - mu / r;
    a = -mu / (2 * energy);
    
    % Результат
    EO = a;
end
