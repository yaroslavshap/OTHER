clear
close all
clc
% ����� ��������� �������
Init = 0;
% �������� ������������
a = 0;
b = 1;
% �������� ��������
c = 1;
% �����
N_grid = 8;
% ������
Km = 0.5;
% ������
for j = 1:N_grid
    % ����� �����
    n(j) = 100 * 2^(j - 1);
    N = n(j);
    % ��� �� ������������
    h = (b - a)/ N;
    % ���������� �����
    x = a + h/2:h:b - h/2;
    % �����
    t = 0;
    t_end = 1;
    % ��� �� ������� � ������ ��� �������
    tau = abs((Km * h) / c);
    K = ((c * tau) / h);
    % ��������� �������
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
    
    % ���������� ���������� �����
    while(t < t_end)
        
        % ��������� �������
        if(c > 0)
            
            %  ������������� ��������� �������
            u_n1(1) = u0(1) - K * (u0(1) - u0(N));
            
            for i = 2:N
                % ����� ����� ������
                u_n1(i) = u0(i) - K * (u0(i) - u0(i - 1));
            end
        end
        
        if(c < 0)
            
            %  ������������� ��������� �������
            u_n1(N) = u0(N) - K * (u0(1) - u0(N));
            
            for i = 1:N - 1
                % ����� ������ ������
                u_n1(i) = u0(i) - K * (u0(i + 1) - u0(i));
            end
        end
        
        % ������� �� ��������� ���� �� �������
        for i = 1:N
            u0(i) = u_n1(i);
        end
        
        t = t + tau;
    end
    
    % ������ �������
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
    
    % �����
    norma_L1(j) = h * sum(abs(u_exact - u_n1));
    
end

% �������� ����������
for j = 2:N_grid
    order(j) = log2(norma_L1(j - 1) / norma_L1(j));
end

% ���������� �������
figure()
plot(x, u_n1, 'b', 'linewidth', 2)
xlabel('x')
ylabel('u(x)')
title('�������')
hold on

plot(x, u_exact, 'k', 'linewidth', 2)
xlabel('x')
ylabel('u(x)')
legend('��������� �������', '������ �������')
grid on

figure()
loglog(n, norma_L1, 'b-o','MarkerSize', 2, 'linewidth', 2)
xlabel('N')
ylabel('||\delta_L_1||')
title('��������� ������ ����������� � �����')
grid on

figure()
semilogy(n, order, 'b-o','MarkerSize', 2, 'linewidth', 2)
xlabel('N')
ylabel('order')
title('��������� ������ �������� ����������')
grid on

