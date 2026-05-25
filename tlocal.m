function t = tlocal(z, TC) % Локальная температура
    dz = z - 125.0;
    if dz <= 0.0
        t = ((-9.8204695e-6 * dz - 7.3039742e-4) * dz^2 + 1.0) * dz * TC(2) + TC(1);
    else
        t = TC(1) + TC(3) * atan(TC(4) * dz * (1 + 4.5e-6 * dz^2.5));
    end
end