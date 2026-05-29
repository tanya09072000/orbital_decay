function [Cx_val, RHO] = Cx()

HEGEO = input('Введите геодезическую высоту, км: ');
Angl_V = input('Угол между нормалью грани и орбитальной скорости, град: ');
userDateStr = input('Введите дату в формате ГГГГ-ММ-ДД: ', 's'); % 2025-01-01

    [Cx_val, RHO, Tinf] = compute_cx( ...
        HEGEO, ...
        Angl_V, ...
        userDateStr);

fprintf('Kоэффициент аэродинамического сопротивления: %.6e\n', Cx_val); % 3.24363433138258
fprintf('Tinf: %.6e\n', Tinf);
end