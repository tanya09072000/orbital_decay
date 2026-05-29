function [Cx_val, RHO, Tinf] = compute_cx( ...
    height_km, angle_deg, date_str, data)
    mu = 398600.4415; % Гравитационный параметр Земли, км³/с²
    Re = 6371; % радиус Земли, км
    aem = 1.6605390666d-27; % атомная единица массы, кг
    k = 1.380649e-23; % постоянная Больцмана, Дж/К

    [TEMP, AMMW, RHO] = cira_main( ...
    height_km, ...
    date_str, ...
    data);
    m = AMMW * aem; % Mолекулярная масса частиц газа (2.65881737978781e-26)
    V = sqrt(mu/(Re+height_km))*1000; % Орбитальная скорость, км/с (7616.56080339665)
    alpha = deg2rad(angle_deg);
    %alpha = Angl_V;
    Tinf = TEMP(1); % 1359.40454756265
    Tw = 0.000201*height_km^3-0.1164*height_km^2+24.217*height_km - 1525; % Температура во фронте ударной волны 6608.5
    S = V/(sqrt((2*k*Tinf)/m)); % 6.4102087037962
    Sw = V/(sqrt((2*k*Tw)/m)); % 2.90733503335382
    S1 = exp(-(S*cos(alpha))^2); % 1.42725932255069e-18
    S2 = sqrt(pi)*S*(cos(alpha))*(1+(1/(2*S^2)))*erf(S*cos(alpha)); % 11.5000515314256
    S3 = ((pi*S)/Sw)*(cos(alpha))^2; %  6.92670928557985
    Cx_val = (2/(sqrt(pi)*S))*(S1+S2+S3);
end

