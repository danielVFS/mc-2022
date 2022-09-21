% Daniel Vitor
% Exercicio 03  

function lista03()
  clc
  format short
  
  A = [-1, 1, -3; 3, -2, 8; 2, -2, 5]; % matriz de coeficiente
  b = [-4; 14; 7]; % vetor da direita
  
  [A, b] = PivotamentoParcial(A, b, size(A));
  Ab = [A b]; % matrix aumentada
  n = size(Ab, 1)
  
  [result, i] = GaussJordan(Ab, n);
  printf("Foram necessarias %d iterações para encontrar x\n", i);
  disp('Resultado = '); disp(result);
end

function [Ab, i] = GaussJordan(Ab, n)
  % colunas
  for i = 1:n
    Ab(i, :) = Ab(i, :) ./ Ab(i, i)% toda a linha divida pela elemento pivô -> fazendo com que o elemento pivô seja 1;
    % linhas
    for j = 1:n
      if j ~= i % não fazer calculos nas diagonais, afinal vao ser 1
        % zerar todos elementos que nao sao as diagonais
        multiplier = Ab(j, i);
        Ab(j, :) = Ab(j, :) - multiplier * Ab(i, :);
      end
    end
  end
end 


function [A, b] = PivotamentoParcial(A, b, n)
  X = zeros(n, 1);

  % funcao de pivotamento
  for k=1:n-1   % n-1 -> na ultima coluna não tem nada pra eliminar abaixo da coluna principal
    max = abs(A(k, k)); % valor da diagonal principal
    p = k; 
    % comparar a diagonal principal com o valor abaixo dela
    for i = k+p:n
      if A(i,k) > max
        max = A(i,k);
        p = i;
      end
    end
    
    if p ~= k % maior elemento não está na diagonal principal
      for j = k:n
        aux =  A(p,j);
        A(p,j) = A(k,j);
        A(k,j) = aux;
      end
      
      aux = b(p);
      b(p) = b(k);
      b(k) = aux;  
    end
  end
end
