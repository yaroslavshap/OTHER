clear
close all
clc
% Выбор начальных условий
Init = 0;
% Интервал пространства
a = 0;
b = 1;
% Скорость переноса
c = 1;
% Сетки
N_grid = 8;
% Курант
Km = 0.5;
% Расчет
for j = 1:N_grid
    % Число ячеек
    n(j) = 100 * 2^(j - 1);
    N = n(j);
    % Шаг по пространству
    h = (b - a)/ N;
    % Координаты ячеек
    x = a + h/2:h:b - h/2;
    % Время
    t = 0;
    t_end = 1;
    % Шаг по времени и курант для расчета
    tau = abs((Km * h) / c);
    K = ((c * tau) / h);
    % Начальные условия
    for i = 1:N
        if(Init == 0)
            u0(i) = 3 * exp(-(20 * x(i) - 8)^2) + 2 * exp(-(26 * x(i) - 13)^2);
        end
        if(Init == 1)
            if ((x(i) < 0.3) || (x(i) > 0.5))
                u0(i) = 0;
            else
                u0(i) = 1;
            end
        end
    end
    
    % Реализация разностной схемы
    while(t < t_end)
        
        % Численное решение
        if(c > 0)
            
            %  Периодические граничные условия
            u_n1(1) = u0(1) - K * (u0(1) - u0(N));
            
            for i = 2:N
                % Явный левый уголок
                u_n1(i) = u0(i) - K * (u0(i) - u0(i - 1));
            end
        end
        
        if(c < 0)
            
            %  Периодические граничные условия
            u_n1(N) = u0(N) - K * (u0(1) - u0(N));
            
            for i = 1:N - 1
                % Явный правый уголок
                u_n1(i) = u0(i) - K * (u0(i + 1) - u0(i));
            end
        end
        
        % Переход на следующий слой по времени
        for i = 1:N
            u0(i) = u_n1(i);
        end
        
        t = t + tau;
    end
    
    % Точное решение
    for i = 1:N
        if(Init == 0)
            u_exact(i) = 3 * exp(-(20 * (mod(x(i) - c * t,b - a) + a) - 8)^2) + 2 * exp(-(26 * (mod(x(i) - c * t,b - a) + a) - 13)^2);
        end
        if(Init == 1)
            if ((mod(x(i) - c * t,b - a) + a < 0.3) || mod(x(i) - c * t,b - a) + a > 0.5)
                u_exact(i) = 0;
            else
                u_exact(i) = 1;
            end
        end
    end
    
    % Норма
    norma_L1(j) = h * sum(abs(u_exact - u_n1));
    
end

% Скорость сходимости
for j = 2:N_grid
    order(j) = log2(norma_L1(j - 1) / norma_L1(j));
end

% Построение решения
figure()
plot(x, u_n1, 'b', 'linewidth', 2)
xlabel('x')
ylabel('u(x)')
title('Решение')
hold on

plot(x, u_exact, 'k', 'linewidth', 2)
xlabel('x')
ylabel('u(x)')
legend('Численное решение', 'Точное решение')
grid on

figure()
loglog(n, norma_L1, 'b-o','MarkerSize', 2, 'linewidth', 2)
xlabel('N')
ylabel('||\delta_L_1||')
title('Численная оценка погрешности в норме')
grid on

figure()
semilogy(n, order, 'b-o','MarkerSize', 2, 'linewidth', 2)
xlabel('N')
ylabel('order')
title('Численная оценка скорости сходимости')
grid on

