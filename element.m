function EO = element(Y1)
% Преобразует вектор состояния в классические кеплеровы элементы орбиты.
% Вход:  Y1 — [x, y, z, Vx, Vy, Vz], км и км/с
% Выход: EO — [a (км), e, i (рад), Omega (рад), omega (рад), nu (рад)]

    C  = orb_constants();
    mu = C.mu;

    r_vec = Y1(1:3)';
    v_vec = Y1(4:6)';
    r     = norm(r_vec);
    v     = norm(v_vec);

    energy = v^2 / 2 - mu / r;
    if energy >= 0
        error('element: незамкнутая орбита (energy = %.4g км²/с²)', energy);
    end
    a = -mu / (2 * energy);

    h_vec = cross(r_vec, v_vec);
    h     = norm(h_vec);

    e_vec = cross(v_vec, h_vec) / mu - r_vec / r;
    e     = norm(e_vec);

    inc = acos(max(-1, min(1, h_vec(3) / h)));

    N_vec = cross([0; 0; 1], h_vec);
    N     = norm(N_vec);
    if N < eps
        Omega = 0;
    else
        Omega = acos(max(-1, min(1, N_vec(1) / N)));
        if N_vec(2) < 0, Omega = 2*pi - Omega; end
    end

    if N < eps || e < eps
        omega = 0;
    else
        omega = acos(max(-1, min(1, dot(N_vec, e_vec) / (N * e))));
        if e_vec(3) < 0, omega = 2*pi - omega; end
    end

    nu = acos(max(-1, min(1, dot(e_vec, r_vec) / (e * r))));
    if dot(r_vec, v_vec) < 0, nu = 2*pi - nu; end

    EO = [a, e, inc, Omega, omega, nu];
end
