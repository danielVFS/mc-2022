% Daniel Vitor
% exercicio 08
function lista04()
  clc;
  format long;
 
  x0 = [3; 0];
  maxIterations = 1000;
  tolerance = 1e-3;
  flag = 1;
  pontoOtimo = [1; 1];

  [xr1, xTodasIteracoes, fxTodasIteracoes, i] = otimizacao(x0, tolerance, maxIterations, pontoOtimo, flag);
  
  printf("O resultado de F(%.6f, %.6f) = %.6f, e foram necessarias %d iteracoes\n", xr1(1), xr1(2), funcao(xr1), i-1);
  
  % -------------------------------------------------------------------------- %

  qtdIteracoes = i-1;
  maiorValor = max(max(xTodasIteracoes));
  menorValor = min(min(xTodasIteracoes));
  intervalo = (menorValor - 1):0.5:(maiorValor + 1);
 
  Z = calculaGrafico3D(intervalo, intervalo);
  
  plotaGrafico(Z, xTodasIteracoes, fxTodasIteracoes ,intervalo, qtdIteracoes);
 
  plotarConvergencia(xTodasIteracoes, fxTodasIteracoes, i-1);
endfunction

function [x0, xTodasIteracoes, fxTodasIteracoes, i] = otimizacao(x0, tolerance, maxIterations, pontoOtimo, flag)
  xTodasIteracoes = zeros(maxIterations, 2); 
  fxTodasIteracoes = zeros(maxIterations, 1);
  
   if flag == 1
     alpha = 0.1;
   else
     alpha = 0.12;
   endif
   
  xTodasIteracoes(1, :) = x0;
  fxTodasIteracoes(1, :) = funcao(x0);
  
  for i = 1:maxIterations
    if flag == 1
     x1 = x0 - alpha*gradiente(x0);
    else
     x1 = x0 - alpha*(hessian_levy13(x0)\gradiente(x0));
    endif
    
    if max(abs(x1 - x0)) < tolerance
      break;
    endif
    
    x0 = x1;
    
    xTodasIteracoes(i, :) = x0;
    fxTodasIteracoes(i, :) = funcao(x0); 
  endfor
  
  distancia_euclidiana = calcula_funcao_euclidiana(x0, pontoOtimo)
  
  xTodasIteracoes = xTodasIteracoes(1:i, :);
  fxTodasIteracoes = fxTodasIteracoes(1:i, :);
endfunction

function d = calcula_funcao_euclidiana(x0, pontoOtimo)
  x1 = x0(1);
  x2 = x0(2);
  y1 = pontoOtimo(1);
  y2 = pontoOtimo(2);
  
  d = sqrt((x1 - y1)^2 + (x2 - y2)^2);
endfunction

function y = funcao(x)
   x1 = x(1);
   x2 = x(2);
   
   y = sin(3*pi*x(1))^2 + (x(1) - 1)^2 * (1 + sin(3*pi*x(2))^2) + (x(2) - 1)^2 * (1 + sin(2*pi*x(2))^2);
endfunction

function g = gradiente(x)
   x1 = x(1);
   x2 = x(2);
    
   g = zeros(2, 1);
   g(1) = 2*sin(3*pi*x(1))*3*pi*(x(1) - 1) + 2*(x(1) - 1);
   g(2) = 2*sin(3*pi*x(2))*3*pi*(x(1) - 1) + 2*(x(2) - 1)*(1 + sin(2*pi*x(2))^2) + 2*(x(2) - 1);
endfunction
 
function h = hessiana(x)
   x1 = x(1);
   x2 = x(2);
   
   h = zeros(2, 2);
   h(1, 1) = 18*pi^2*sin(3*pi*x(1))^2 + 2*(1 + sin(3*pi*x(2))^2);
   h(1, 2) = 6*pi*sin(3*pi*x(2))*sin(3*pi*x(1));
   h(2, 1) = 6*pi*sin(3*pi*x(2))*sin(3*pi*x(1));
   h(2, 2) = 4*(x(2) - 1)*sin(2*pi*x(2))^2 + 2*(1 + sin(2*pi*x(2))^2);
endfunction

function plotarConvergencia(xTodasIteracoes, fxTodasIteracoes, i)
  figure(2);
  clf;
  
  subplot(2, 1, 1);
  p3 = plot(0:(i-1), xTodasIteracoes(1:i, :), 'linewidth', 2);
  xlim([0, i-1]);
  set(gca, 'fontsize', 20);
  grid('on');
  title('Grafico de convergencia de x1 e x2');
  xlabel('Iteracoes');
  ylabel('x');
  
  # -------------------------------------------------------------------------- #

  subplot(2, 1, 2);
  p4 = plot(0:(i-1), fxTodasIteracoes(1:i, :), 'linewidth', 2);
  xlim([0, i-1]);
  set(gca, 'fontsize', 20);
  grid('on');
  title("Convergencia de F(xr)");
  xlabel("Iteracoes");
  ylabel("F(x)");
endfunction

function plotaGrafico(Z, X, Y ,intervalo, qtdIteracoes)
  figure(1);
  clf;
  
   for i = 1:qtdIteracoes-1
    pause(0);
    s1 = surf(intervalo, intervalo, Z(:, :, 1));
    view([-10 -90 90]);
    hold on
    plot3(X(1:i, 2), X(1:i, 1), Y(1:i, 1), 'color', 'm', 'linewidth', 2);
    p1 = plot3(X(i, 2), X(i, 1), Y(i, 1), 'Marker', 'o', 'MarkerFaceColor', 'm', 'color', 'k', 'linewidth', 2, 'MarkerSize', 10);
    set(gca, 'fontsize', 15);
    hold off
    title(sprintf('f(%.6f, %.6f) = %.6f, Iteracoes: %i', X(i, 1), X(i, 2), Y(i, 1), i));
    xlabel('y');
    ylabel('x');
    zlabel('u(x1, x2)');
    legend([s1, p1], {'U(x1, x2)', 'Posicao calculada'}, 'location', 'northeast');
  endfor
  
end

function Z = calculaGrafico3D(intervaloX, intervaloY)
  numX = numel(intervaloX);
  numY = numel(intervaloY);
  
  Z = zeros(numX, numY, 1);
  
  for i = 1:numX
    for j = 1:numY
      Z(i, j, 1) = funcao([intervaloX(i); intervaloY(j)]);
    endfor
   endfor
end
  
