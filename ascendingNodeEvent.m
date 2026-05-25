function [value, isterminal, direction] = ascendingNodeEvent(~, y)
    value = y(3);         % z-компонента = 0 — экватор
    isterminal = 1;       % остановить интегрирование
    direction = +1;       % пересекаем снизу вверх (восходящий узел)
end