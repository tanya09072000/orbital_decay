function dy_dt = two_body(t, y, C, Sm, m_total, h_grid, rho_grid, Cx_grid)
% Правая часть уравнения движения: гравитация + J2 + аэродинамическое торможение.
% Вход:  y        — [x, y, z, Vx, Vy, Vz], км и км/с
%        C        — структура физических констант (orb_constants)
%        Sm       — площадь миделевого сечения, м²
%        m_total  — масса КА, кг
%        h_grid   — сетка высот для интерполяции, км
%        rho_grid — плотность атмосферы на сетке h_grid, кг/м³
%        Cx_grid  — коэффициент сопротивления на сетке h_grid

    mu    = C.mu;
    J2    = C.J2;
    Re    = C.Re;
    omega = C.omega;

    r_vec  = y(1:3);
    v_vec  = y(4:6);
    r_norm = norm(r_vec);

    r_si = r_vec * 1000;   % позиция в м
    v_si = v_vec * 1000;   % скорость в м/с

    zr_ratio = (r_vec(3) / r_norm)^2;

    v_atm      = cross([0; 0; omega], r_si);   % скорость атмосферы, м/с
    v_rel      = v_si - v_atm;                 % относительная скорость, м/с
    v_rel_norm = norm(v_rel);

    a_grav = (-mu / r_norm^3) * r_vec;

    h_km   = r_norm - Re;
    rho    = interp1(h_grid, rho_grid, h_km, 'linear', 0);
    Cx_val = interp1(h_grid, Cx_grid,  h_km, 'linear', Cx_grid(end));

    a_drag = (-0.5 * Cx_val * Sm * rho * v_rel_norm / m_total / 1000) * v_rel;

    a_J2 = ((-3 * J2 * mu * Re^2) / (2 * r_norm^5)) * ...
           [ r_vec(1) * (1 - 5*zr_ratio);
             r_vec(2) * (1 - 5*zr_ratio);
             r_vec(3) * (3 - 5*zr_ratio) ];

    dy_dt = [v_vec; a_grav + a_drag + a_J2];
end
