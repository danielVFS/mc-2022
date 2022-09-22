% Daniel Vitor
% Exercicio 11  

function lista03()
  clc
  format short
 
  A = [-1, 1, -3; 3, -2, 8; 2, -2, 5]; % matriz de coeficiente
  a = [-4; 14; 7]; % vetor da direita
   
  Aa = PivotamentoParcial([A a], size(A)); % [A b] -> matrix aumentada
  n = size(Aa, 1)
 
  [result, i] = GaussJordan(Aa, n);
  printf("Foram necessarias %d iterações para encontrar x\n", i);
  disp('Resultado = '); disp(result);
end

function [Aa, i] = GaussJordan(Aa, n)
  % colunas
  for i = 1:n
    Aa(i, :) = Aa(i, :) ./ Aa(i, i);% toda a linha divida pela elemento pivô -> fazendo com que o elemento pivô seja 1;
    % linhas
    for j = 1:n
      if j ~= i % não fazer calculos nas diagonais, afinal vao ser 1
        % zerar todos elementos que nao sao as diagonais
        multiplier = Aa(j, i);
        Aa(j, :) = Aa(j, :) - multiplier * Aa(i, :);
      end
    end
  end
end


function Aa = PivotamentoParcial(Aa, n)
  for i = 1:n-1
      [~, pos] = max(abs(Aa(i:n, i))); % Ab(col:end, col)
      pos = pos + i -1;
      if pos ~= i
        aux = Aa(pos, :);
        Aa(pos, :) = Aa(i, :);
        Aa(i, :) = aux;
      end
  end
 end