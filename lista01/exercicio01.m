# Daniel Vitor
# exercicio - 09

function exercicio01()
  # convenção de variaveis seguindo o apresentado em aula/slides.
  xu = 1; # interval1
  xl = 0; # interval2
   
  tolerance = 1e-5;
  
  xr0 = inf;

  vetx = zeros(1,1000);
  vety = zeros(1,1000);
  
  x = -2:0.05:3;
  y = zeros(size(x));
  
  counter = 1;
  for getValues = x
    y(counter) = getRoot(getValues);
    counter++;
  endfor
  
  maxIterations = 1000;
  
  for iterations = 1:maxIterations
    # Fórmula que calcula a posição da raiz
    xr = xu - ((getRoot(xu)*(xl-xu))/(getRoot(xl)-getRoot(xu)));
    
    if getRoot(xl)*getRoot(xr) > 0
      xl = xr; # Raíz na primeira metade do intervalo 
    else
      xu = xr; # Raíz na segunda metade do intervalo
    endif   
    
    plotGraph(x, y, iterations, xr);
    
    if abs(xr0 - xr) <= tolerance
      break;
    endif
    
    xr0 = xr;
    
    vetx(iterations) = xr0;
    vety(iterations) = getRoot(xr0);
  end
  
  printf("F(%.6f) = %.6f, sendo necessario %d iteracoes\n", xr, getRoot(xr), iterations);
  vetx = vetx(1:iterations);
  vety = vety(1:iterations);
  
  plotConvergence(vetx, vety, iterations);
end

function y = getRoot(x)
  y = (x.^3)+(2.*(x.^2))-2;
end

function plotGraph(x, y, i, xr)
    figure(1)
    clf
    grid('on');
    hold('on');
    
    p1 = plot(x, y, 'linewidth', 2, 'og');
    p2 = plot(xr, getRoot(xr), 'linewidth', 2, 'oxr', 'markersize', 15);
    
    hold('off');
    
    # Informações do Gráfico
    
    set(gca, 'fontsize', 15);
    title(sprintf('F(%.6f) = %.6f, Iteracoes: %d', xr, getRoot(xr), i));
    xlabel('x');
    ylabel('F(x)');
    legend([p1, p2], {'Funcao', 'Raiz calculada'}, 'location', 'northeast');
    pause(0.10); # Velocidade das iteracoes
end

function plotConvergence(vetx, vety, i)
  figure(2)
  clf
  
  subplot(2, 1, 1);
  p3 = plot(1:(i-1), vetx(1:(i-1)), 'linewidth', 2, 'g');
  set(gca, 'fontsize', 20);
  grid('on');
  title("Convergencia de Xr");
  xlabel("Iteracao");
  ylabel("x");
  
  subplot(2, 1, 2);
  p4 = plot(1:(i-1), vety(1:(i-1)), 'linewidth', 2, 'color', 'b');
  set(gca, 'fontsize', 20);
  grid('on');
  title("Convergencia de F(xr)");
  xlabel("Iteracao");
  ylabel("F(x)");
end