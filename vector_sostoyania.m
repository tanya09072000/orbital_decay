function VS = vector_sostoyania(input_data)
% Переводит кеплеровы элементы орбиты в вектор состояния.
% Вход:  [a (км), e, i (°), Omega (°), omega (°), nu (°)]
% Выход: [x, y, z (км), Vx, Vy, Vz (км/с)]

    C  = orb_constants();
    mu = C.mu;

    a   = input_data(1);
    e   = input_data(2);
    i   = deg2rad(input_data(3));
    Om  = deg2rad(input_data(4));
    w   = deg2rad(input_data(5));
    nu  = deg2rad(input_data(6));

    uu = w + nu;
    p  = a * (1 - e^2);
    r  = p / (1 + e * cos(nu));

    Vr = sqrt(mu / p) * e * sin(nu);
    Vu = sqrt(mu * p) / r;

    x  =  r * (cos(Om)*cos(uu) - sin(Om)*sin(uu)*cos(i));
    y  =  r * (sin(Om)*cos(uu) + cos(Om)*sin(uu)*cos(i));
    z  =  r *  sin(uu) * sin(i);

    Vx = Vr*(cos(Om)*cos(uu) - sin(Om)*sin(uu)*cos(i)) - Vu*(cos(Om)*sin(uu) + sin(Om)*cos(uu)*cos(i));
    Vy = Vr*(sin(Om)*cos(uu) + cos(Om)*sin(uu)*cos(i)) - Vu*(sin(Om)*sin(uu) - cos(Om)*cos(uu)*cos(i));
    Vz = Vr*sin(uu)*sin(i) + Vu*cos(uu)*sin(i);

    VS = [x, y, z, Vx, Vy, Vz];
end
