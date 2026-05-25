function VS = vector_sostoyania(input_data) 
% перевод элементов орбиты в вектор состояния КА
% ВХОД: большая полуось (км), эксцентриситет, наклонение (град), долгота восходящего
% узла (град), аргумент перицентра (град), истинная аномалия (град)
% гравитационный параметр (км^3/с^2)
% ВЫХОД: вектор положения (км), вектор скорости (км/с) 

a = input_data (1);
e = input_data (2);
i = input_data (3);
om = input_data (4);
w = input_data (5);
u = input_data (6);
mu = double(398600.4415);

% преобразование углов в радианы
i = deg2rad(i);
om = deg2rad(om);
w = deg2rad(w);
u = deg2rad(u);

uu=w+u;

% расчет расстояния от объекта до центрального тела
p = a*(1-e^2);
r = p/(1+e*cos(u));

Vr = sqrt(mu/p)*e*sin(u);
Vu = sqrt((mu*p)/(r*r));

% компоненты вектора состояния
x = r * (cos(om)*cos(uu)-sin(om)*sin(uu)*cos(i));
y = r * (sin(om)*cos(uu)+cos(om)*sin (uu)*cos(i));
z = r * (sin(uu)*sin(i));

% компоненты вектора скорости
Vx = Vr*(cos(om)*cos(uu)-sin(om)*sin(uu)*cos(i))-Vu*(cos(om)*sin(uu)+sin(om)*cos(uu)*cos(i));
Vy = Vr*(sin(om)*cos(uu)+cos(om)*sin(uu)*cos(i))-Vu*(sin(om)*sin(uu)-cos(om)*cos(uu)*cos(i));
Vz = Vr*sin(uu)*sin(i)+Vu*cos(uu)*sin(i);

VS = [x, y, z, Vx, Vy, Vz];
end