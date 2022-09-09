% Daniel Vitor
% Exercicio 03

function lista02()
  clc
  format long
  
  x0 = [1; 1];
  
  tol = 1e-5;
  
  X = zeros(1001, 2); % x,y
  Y = zeros(1001, 2); % u(x,y), v(x,y)
  
  X(1, 1:2) = x0;
  Y(1, 1:2) = sistema(x0);
  
  
  for i = 2:1001
    J = jacobiano(x0);
    f = sistema(x0);
  
    x1 = x0 - inv(J)*f;
    
    X(i+1, 1:2) = x1;
    Y(i+1, 1:2) = sistema(x1);
           
    if max(abs(x1 - x0)) <= tol
      x0 = x1;
      break;
    endif
    
    x0 = x1;
  end
  
  value = sistema(x1);
  printf("O resultado de F(%.6f) = %.6f, e foram necessorias %d iterações\n", x1(1), value(1), i-1);
  printf("O resultado de F(%.6f) = %.6f, e foram necessorias %d iterações\n", x1(2), value(2), i-1);
  
  % only the ones i want
  X = X(1:(i+1), :);
  Y = Y(1:(i+1), :);
  
  qtdIteracoes = i;
  maiorValor = max(max(X));
  menorValor = min(min(X));
  intervalo = (menorValor - 1):0.1:(maiorValor + 1);
  
  Z = calculaGrafico3D(intervalo, intervalo);
  
  plotaGrafico(Z, X, Y ,intervalo, qtdIteracoes);
  
  plotarConvergencia(X, Y, i);
end

function plotaGrafico(Z, X, Y ,intervalo, qtdIteracoes)
  figure(1);
  clf;
  
  for i = 1:qtdIteracoes
    pause(1);
    subplot(2,1,1);
    s1 = surf(intervalo, intervalo, Z(:, :, 1));
    view([-10 -90 90]);
    hold on
    plot3(X(1:i, 2), X(1:i, 1), Y(1:i, 1), 'color', 'm', 'linewidth', 2);
    p1 = plot3(X(i, 2), X(i, 1), Y(i, 1), 'Marker', 'o', 'MarkerFaceColor', 'm', 'color', 'k', 'linewidth', 2, 'MarkerSize', 10);
    set(gca, 'fontsize', 15);
    hold off
    title(sprintf('f(%.6f) = %.6f, Iteracoes: %i', X(i, 1), Y(i, 1), i-1));
    xlabel('y');
    ylabel('x');
    zlabel('u(x1, x2)');
    legend([s1, p1], {'U(x1, x2)', 'Posicao calculada'}, 'location', 'northeast');
        
    subplot(2,1,2);
    s2 = surf(intervalo, intervalo, Z(:, :, 2));
    view([-10 -90 90]);
    hold on
    plot3(X(1:i, 2), X(1:i, 1), Y(1:i, 2), 'color', 'm', 'linewidth', 2);
    p2 = plot3(X(i, 2), X(i, 1), Y(i, 2), 'Marker', 'o', 'MarkerFaceColor', 'm', 'color', 'k', 'linewidth', 2, 'MarkerSize', 10);
    set(gca, 'fontsize', 15);
    hold off
    title(sprintf('f(%.6f) = %.6f, Iteracoes: %i', X(i, 2), Y(i, 2), i-1));
    xlabel('y');
    ylabel('x');
    zlabel('v(x1, x2)');
    legend([s2, p2], {'V(x1, x2)', 'Posicao calculada'}, 'location', 'northeast');
  endfor
end

function Z = calculaGrafico3D(intervaloX, intervaloY)
  numX = numel(intervaloX);
  numY = numel(intervaloY);
  
  Z = zeros(numX, numY, 2);
  
  for i = 1:numX
    for j = 1:numY
      Z(i, j, 1:2) = sistema([intervaloX(i); intervaloY(j)]);
    endfor
   endfor
end

function vetY = sistema(vetX)
  x1 = vetX(1);
  x2 = vetX(2);
  
  vetY = zeros(2,1);
  vetY(1) = 2.*x1 - 4.*x1.*x2 + 2.*(x2.^2);
  vetY(2) = 3.*(x2.^2) + 6.*x1 - x1.^2 - 4.*x1.*x2 - 5;
end

function J = jacobiano(vetX)
  x1 = vetX(1);
  x2 = vetX(2);
  
  J = zeros(2,2);
  J(1,1) = -4.*x2 + 2;         % d u(x1, x2)/x1
  J(1,2) = 4.*x2 - 4.*x1;      % d u(x1, x2)/x2
  J(2,1) = 6 - 4.*x2 - 2.*x1;  % d v(x1, x2)/x1
  J(2,2) = 6.*x2 - 4.*x1;      % d v(x1, x2)/x2
end

function plotarConvergencia(X, Y, i)  
  figure(2)
  clf
  subplot(2, 1, 1);
  p3 = plot(0:(i-1), X(1:(i), :), 'linewidth', 2);
  
  set(gca, 'fontsize', 20);
  grid('on');
  title("Gráfico de convergência de x1 e x2");
  xlabel("Iterações");
  ylabel("x"); 
  
  # -------------------------------------------------------------------------- #
  
  subplot(2, 1, 2);
  p4 = plot(0:(i-1), Y(1:(i), :), 'linewidth', 2);
  set(gca, 'fontsize', 20);
  grid('on');
  title("Grafico de convergencia de F(x1 e x2)");
  xlabel("Iterações");
  ylabel("F(x)");
end

