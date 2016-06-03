
stacksize('max');
clear();
funcprot(0);
printf("\n\n--- (©) Alexei Dergunov, 2010 ---\n");

//----------------------------------------------
printf("\n----------- Задание 1 -----------\n");
//----------------------------------------------

x = fscanfMat('C:\data2.txt'); // считываем выборку из файла

n = length(x); // длина выборки
if n <> 5000 then
  printf("Выборка неверно считана из файла!\n");
end;

// считает корреляционную функцию
function[ans] = corr(x, t)
  ans = 0;
  k = abs(t); // т.к. корреляционная функция четная, R(-t) = R(t)
  n = length(x);
  mx = mean(x);
  // просто считаем по формуле
  for i = 1 : n - k
    ans = ans + (x(i) - mx) * (x(i + k) - mx);
  end;
  ans = ans / (n - k);
endfunction

// считает нормированную корреляционную функцию
function[ans] = normCorr(x, t)
  n = length(x);
  dx = corr(x, 0);
  // делим на дисперсию
  ans = corr(x, t) / dx;
endfunction

// считает значения корреляционной функции
// R(t), t = [0, tmax]
function[ans] = getCorrValues(x, tmax)
  ans = zeros(tmax + 1, 1);
  for i = 0 : tmax
    ans(i + 1) = corr(x, i);
  end;
endfunction

// считает значения нормированной корреляционной функции
// r(t), t = [0, tmax]
function[ans] = getNormCorrValues(x, tmax)
  ans = zeros(tmax + 1, 1);
  for i = 0 : tmax
    ans(i + 1) = normCorr(x, i);
  end;
endfunction

// считает радиус корреляции
function[ans] = corrRadius(x, tmax)
  // начинаем проверять радиус корреляции
  // с t = tmax в меньшую сторону
  for t = tmax : -1 : 0
    if abs(normCorr(x, t)) >= 1 / %e
      // как только появилось значение, большее 1 / e,
      // выходим из цикла
      ans = t + 1;
      break;
    end;
  end;
endfunction

// строит график нормированной корреляционной функции, t = [0, tmax]
function[] = plotNormCorrFunc(x, tmax, _style)
  t = [0 : tmax]; // значения по оси абсцисс
  r = getNormCorrValues(x, tmax); // значения по оси ординат
  plot2d(t, r, style = _style);
endfunction

mx = mean(x); // выборочное среднее
printf("Выборочное среднее   = %.4f\n", mx);

dx = corr(x, 0); // выборочная дисперсия
printf("Выборочная дисперсия = %.4f\n", dx);

T = corrRadius(x, 250); // радиус корреляции
printf("Радиус корреляции    = %d\n", T);


//----------------------------------------------
printf("\n----------- Задание 2 -----------\n");
//----------------------------------------------

// находит параметры beta(1)...beta(M) модели АРСС
function[ans] = findBetas(x, M, N)
  if M == 0 then
    ans = [];
  else
    R = zeros(2 * M, 1);
    // далее будут хитрые, но верные преобразования индексов
    for i = N - M + 1 : N + M
      // нужные нам значения корреляционной функции
      R(i + M - N) = corr(x, i);
    end;
    A = zeros(M, M); // матрица системы
    b = zeros(M, 1); // столбец свободных членов
    for i = 1 : M
      b(i) = -R(M + i);
      for j = 1 : M
        A(i, j) = R(M + i - j);
      end;
    end;
    ans = linsolve(A, b); // решаем систему Ax + b = 0
  end;
endfunction

// находит параметры alpha(0)...alpha(N) модели АРСС
// по методу 2
function[ans] = findAlphas2(x, M, N, betas)
  
  // вспомогательная функция для метода 2
  // надо сгенерировать модель СС(N)
  // для последовательности y,
  // используя исходную последовательность x
  function[ans] = generateMovingAverage(x, betas)
    n = length(x);
    ans = zeros(n, 1);
    M = length(betas);
    // просто считаем по формуле
    for i = 1 : n
      ans(i) = x(i);
      j = 1;
      while j <= M & i - j > 0 do
        ans(i) = ans(i) - betas(j) * x(i - j);
        j = j + 1;
      end;
    end;
  endfunction
  
  // генерируем последовательность y
  y = generateMovingAverage(x, betas);
  
  // находим значения ее корреляционной функции R(n), n = [0, N]
  R = getCorrValues(y, N);
  
  // функция, необходимая для решения нелинейной системы
  // она должна быть равна нулю ...
  function[ans] = funcToZero(alpha)
    ans = zeros(N + 1, 1);
    for i = 0 : N
      ans(i + 1) = -R(i + 1); // ... поэтому надо написать так
      for j = 0 : N - i
        ans(i + 1) = ans(i + 1) + alpha(j + 1) * alpha(j + i + 1);
      end;
    end;
  endfunction
  
  // решаем нелинейную систему
  [ans, values, info] = fsolve([1 : N + 1], funcToZero);
  if info == 4 then // итерации завершились неудачно - решения не существует
    ans(1) = %nan;
  end;
  
endfunction

// взаимная корреляционная функция
// при использовании метода 1
// эта функция нужна для вычисления alpha
function[ans] = mutCorr(alphas, betas)
  N = length(alphas) - 1;
  M = length(betas);
  ans = zeros(N + 1, 1);
  // считаем по формуле
  for i = 0 : N
    ans(i + 1) = alphas(i + 1);
    m = min(i, M);
    for j = 1 : m
      ans(i + 1) = ans(i + 1) + betas(j) * ans(i - j + 1);
    end;
  end;
endfunction

// находит параметры alpha(0)...alpha(N) модели АРСС
// по методу 1
function[ans] = findAlphas1(x, M, N, betas)
  // находим значения корреляционной функции R(m), m = [0, max(M, N)]
  R = getCorrValues(x, max(M, N));

  // функция, необходимая для решения нелинейной системы
  // она должна быть равна нулю
  function[ans] = funcToZero(alphas)
    ans = zeros(N + 1, 1);
    for i = 0 : N
      ans(i + 1) = -R(i + 1); // "свободный член"
      // слагаемые, содержащие beta
      for j = 1 : M
        ans(i + 1) = ans(i + 1) + betas(j) * R(abs(j - i) + 1);
      end;
      RR = mutCorr(alphas, betas); // взаимная корреляционная функция
      // слагаемые, содержащие alpha
      for j = i : N
        ans(i + 1) = ans(i + 1) + alphas(j + 1) * RR(j - i + 1);
      end;
    end;
  endfunction

  // решаем нелинейную систему
  [ans, values, info] = fsolve([1 : N + 1], funcToZero);
  if info == 4 then // итерации завершились неудачно - решения не существует
    ans(1) = %nan;
  end;
  
endfunction

// проверяет устойчивость модели
function[ans] = isStable(betas)
  n = length(betas);
  select n
    case 0 then
      ans = %T;
    case 1 then
      ans = abs(betas(1)) < 1;
    case 2 then
      ans = abs(betas(2)) < 1 &..
            abs(betas(1)) < 1 - betas(2);
    case 3 then
      ans = abs(betas(3)) < 1 &..
            abs(betas(1) + betas(3)) < 1 - betas(2) &..
            abs(betas(2) + betas(1) * betas(3)) < abs(1 - betas(3) ^ 2);
  end;
endfunction

// красивый форматированный вывод вектора в строчку
function[] = printfVector(str, x)
  n = length(x);
  printf(str);
  for i = 1 : n
    if x(i) >= 0 then
      printf(" ");
    end;
    printf("%.4f", x(i));
    if i < n then
      printf("\t");
    else
      printf("\n");
    end;
  end;
  if n == 0 then
    printf("\n");
  end;
endfunction

function[] = dispAllModels
  for M = 0 : 3
    for N = 0 : 3
      printf("(M, N) = (%d, %d)\n", M, N);
      betas = findBetas(x, M, N);
      if isStable(betas) == %F then
        printf("Модель неустойчива\n");
      else
        alphas = findAlphas1(x, M, N, betas);
        if isnan(alphas(1)) == %T then
          printf("Модели не существует\n");
        else
          printfVector("alphas = ", alphas);
          printfVector("betas  = ", betas);
        end;
      end;
      printf("\n");
    end;
  end;
endfunction

// выводим на экран параметры всех моделей
dispAllModels;


//----------------------------------------------
printf("\n----------- Задание 3 -----------\n");
//----------------------------------------------

// рассчитывает теоретическую корреляционную функцию,
// исходя из системы уравнений Юла - Уокера
function[ans] = theorCorr(alphas, betas, t)
  M = length(betas);
  N = length(alphas) - 1;
  n = max(M + 1, N + 1); // количество уравнений
  RR = mutCorr(alphas, betas);
  b = zeros(n, 1); // столбец свободных членов
  for i = 0 : N
    for j = i : N
      b(i + 1) = b(i + 1) + alphas(j + 1) * RR(j - i + 1);
    end;
  end;
  A = zeros(n, n); // матрица системы
  c = zeros(n, 1); // коэффициенты при R(i) в исходной системе
  c(1) = -1;
  c(2 : M + 1) = betas;
  c(M + 2 : n) = 0;
  // составляем матрицу системы
  for i = 1 : n
    p = i; // индекс, движущийся вправо
    k = i; // индекс, движущийся влево
    for j = 1 : n
      if k > 0 then
      // если в c(i) i > 0
        A(i, j) = A(i, j) + c(k);
      end;
      if p <= n & p <> k then
      // если в c(i) i <= n
      // и индексы, движущиеся влево и вправо, не совпадают
        A(i, j) = A(i, j) + c(p);
      end;
      p = p + 1;
      k = k - 1;
    end;
  end;
  R = linsolve(A, b);
  if t + 1 < n then
    // все запрошенные значения нам уже известны
    ans = R(1 : t + 1);
  else
    // надо вычислить еще немного, по формулам из блока Б
    ans = zeros(t + 1, 1);
    ans(1 : n) = R(1 : n);
    for i = n + 1 : t + 1
      for j = 1 : M
        ans(i) = ans(i) + betas(j) * ans(i - j);
      end;
    end;
  end;
endfunction

function[ans] = theorNormCorr(alphas, betas, t)
  ans = theorCorr(alphas, betas, t);
  dx = ans(1); // дисперсия
  ans = ans / dx; // делим на нее
endfunction

// строит график теоретической нормированной
// корреляционной функции, t = [0, tmax]
function[] = plotTheorNormCorr(alphas, betas, tmax, _style)
  t = [0 : tmax]; // значения по оси абсцисс
  r = theorNormCorr(alphas, betas, tmax); // значения по оси ординат
  plot2d(t, r, style = _style);
endfunction

// среднее квадратическое отклонение
// первых m + 1 отсчетов (т.е. [0, m])
// нормированных корреляционных функций r1 и r2
function[ans] = deviation(r1, r2, m)
  ans = 0;
  for i = 1 : m + 1
    ans = ans + (r1(i) - r2(i)) ^ 2;
  end;
endfunction

// ищет лучшие модели АР, СС и АРСС
function[EpsMatr, bestARparams, bestMAparams, bestARMAparams] = ..
        findBestModels(Mmax, Nmax)

  // находит все отклонения
  EpsMatr = zeros(Mmax + 1, Nmax + 1); // матрица отклонений
  for M = 0 : Mmax
    for N = 0 : Nmax
      betas = findBetas(x, M, N);
      if isStable(betas) == %F then
        EpsMatr(M + 1, N + 1) = %inf; // модель неустойчива
      else
        alphas = findAlphas1(x, M, N, betas);
        if isnan(alphas(1)) == %T then
          EpsMatr(M + 1, N + 1) = %inf; // модели не существует
        else
          // сравниваем теоретическую и выборочную
          // нормированные корреляционные функции
          theor_r = theorNormCorr(alphas, betas, 10);
          r = getNormCorrValues(x, 10);
          EpsMatr(M + 1, N + 1) = deviation(r, theor_r, 10);
        end;
      end;
    end;
  end;

  // находит лучшую модель АР
  bestARparams = zeros(1, 2); // параметры M, N лучшей модели АР  
  AReps = EpsMatr(1, 1);
  for M = 1 : Mmax
    if EpsMatr(M + 1, 1) < AReps then
      AReps = EpsMatr(M + 1, 1);
      bestARparams(1) = M;
    end;
  end;

  // находит лучшую модель СС
  bestMAparams = zeros(1, 2); // параметры M, N лучшей модели СС
  MAeps = EpsMatr(1, 1);
  for N = 1 : Nmax
    if EpsMatr(1, N + 1) < MAeps then
      MAeps = EpsMatr(1, N + 1);
      bestMAparams(2) = N;
    end;
  end;

  // находит лучшую модель АРСС
  bestARMAparams = zeros(1, 2); // параметры M, N лучшей модели АРСС
  if AReps < MAeps then
    ARMAeps = AReps;
    bestARMAparams = bestARparams;
  else
    ARMAeps = MAeps;
    bestARMAparams = bestMAparams;
  end;
  for M = 1 : Mmax
    for N = 1 : Nmax
      if EpsMatr(M + 1, N + 1) < ARMAeps then
        ARMAeps = EpsMatr(M + 1, N + 1);
        bestARMAparams = [M, N];
      end;
    end;
  end;

endfunction

// находим лучшие модели
[EpsMatr, bestARparams, bestMAparams, bestARMAparams] = ..
        findBestModels(3, 3);

// красивый форматированный вывод матрицы
function[] = printfMatrix(str, A)
  [m, n] = size(A);
  printf(str + "\n");
  for i = 1 : m
    for j = 1 : n
      if A(i, j) >= 0 then
        printf(" ");
      end;
      printf("%.4f", A(i, j));
      if j < n then
        printf("\t");
      else
        printf("\n");
      end;
    end;
  end;
endfunction

printfMatrix("Теоретические погрешности моделей АРСС(M, N):", EpsMatr);
printf("\nПараметры наилучшей модели АР (M, N)   = (%d, %d)\n", ..
        bestARparams(1), bestARparams(2));
printf("Параметры наилучшей модели CC (M, N)   = (%d, %d)\n", ..
        bestMAparams(1), bestMAparams(2));
printf("Параметры наилучшей модели АРCC (M, N) = (%d, %d)\n", ..
        bestARMAparams(1), bestARMAparams(2));

// и считаем их параметры
// АР
betasAR = ..
  findBetas(x, bestARparams(1), bestARparams(2));
alphasAR = ..
  findAlphas1(x, bestARparams(1), bestARparams(2), betasAR);
// СС
betasMA = ..
  findBetas(x, bestMAparams(1), bestMAparams(2));
alphasMA = ..
  findAlphas1(x, bestMAparams(1), bestMAparams(2), betasMA);
// АРСС
betasARMA = ..
  findBetas(x, bestARMAparams(1), bestARMAparams(2));
alphasARMA = ..
  findAlphas1(x, bestARMAparams(1), bestARMAparams(2), betasARMA);


//----------------------------------------------
printf("\n----------- Задание 4 -----------\n");
//----------------------------------------------

// спектральная плотность мощности модели АРСС
// с параметрами alphas, betas Fi(exp(iw))
function[ans] = ARMASpecPowDens(alphas, betas, w)
  N = length(alphas) - 1;
  M = length(betas);
  a = 0; // числитель дроби
  for k = 0 : N
    a = a + alphas(k + 1) * exp(%i * w * k);
  end;
  b = 1; // знаменатель дроби
  for k = 1 : M
    b = b - betas(k) * exp(%i * w * k);
  end;
  // еще делим на дисперсию
  ans = (abs((a / b)) ^ 2) / theorCorr(alphas, betas, 0);
endfunction

// оценка спектральной плотности мощности выборки,
// для нахождения используется вектор R отсчетов КФ
function[ans] = specPowDens(R, w)
  n = length(R);
  // просто считаем по формуле
  ans = R(1);
  for k = 1 : n - 1
    ans = ans + 2 * R(k + 1) * cos(w * k);
  end;
  // делим на дисперсию
  ans = ans / R(1);
endfunction

// получает вектор значений спектральной плотности
// мощности модели АРСС с коэффициентами alphas, betas
// Fi(w), w = [0, pi] с шагом step
function[ans] = getARMASpecPowDensValues(alphas, betas, step)
  n = int(%pi / step) + 1;
  ans = zeros(n, 1);
  w = 0;
  for i = 1 : n
    ans(i) = ARMASpecPowDens(alphas, betas, w);
    w = min(w + step, %pi);
  end;
endfunction

// получает вектор значений спектральной плостности
// мощности выборки x Fi(w), w = [0, pi] с шагом step,
// используя отсчеты R(t), t = [0, tmax]
function[ans] = getSpecPowDensValues(x, tmax, step)
  n = int(%pi / step) + 1;
  R = getCorrValues(x, tmax);
  ans = zeros(n, 1);
  w = 0;
  for i = 1 : n
    ans(i) = specPowDens(R, w);
    w = min(w + step, %pi);
  end;
endfunction

// строит график спектральной плостности
// мощности модели АРСС
function[] = plotARMASpecPowDens(alphas, betas, step, _style)
  w = [0 : step : %pi];
  Fi = getARMASpecPowDensValues(alphas, betas, step);
  plot2d(w, Fi, style = _style);
endfunction

// строит график оценки спектральной
// плостности мощности выборки x
function[] = plotSpecPowDens(x, tmax, step, _style)
  w = [0 : step : %pi];
  Fi = getSpecPowDensValues(x, tmax, step);
  plot2d(w, Fi, style = _style);
endfunction


//----------------------------------------------
printf("\n----------- Задание 5 -----------\n");
//----------------------------------------------

// генерирует count значений случайного процесса АРСС
function[ans] = generateARMA(alphas, betas, count)
  badCount = 1000;
  n = count + badCount;
  ans = zeros(n, 1);
  M = length(betas);
  N = length(alphas) - 1;
  ksi = grand(n, 1, 'nor', 0, 1); // ksi(i) ~ N(0, 1)
  // просто считаем по формуле
  for i = 1 : n
    ans(i) = 0;
    // добавляем N + 1 слагаемое с alpha(j)
    j = 0;
    while j <= N & i - j > 0 do
      ans(i) = ans(i) + alphas(j + 1) * ksi(i - j);
      j = j + 1;
    end;
    // добавляем M слагаемых с beta(j)
    j = 1;
    while j <= M & i - j > 0 do
      ans(i) = ans(i) + betas(j) * ans(i - j);
      j = j + 1;
    end;
  end;
  // отбрасываем первые badCount = 1000 отсчетов
  ans = ans(badCount + 1 : n);
  // и прибавляем мат. ожидание mx
  ans = ans + mx;
endfunction

// генерируем 5000 значений каждой модели
imitAR = generateARMA(alphasAR, betasAR, 5000);
imitMA = generateARMA(alphasMA, betasMA, 5000);
imitARMA = generateARMA(alphasARMA, betasARMA, 5000);

// строит реализацию count значений случайного процесса
// если params = true, строятся также линии среднего значения и СКО
function[] = plotRealization(x, count, _style, params)
  n = min(count, length(x));
  // строим саму реализацию
  plot2d([1 : n], x(1 : n), style = _style);
  if params == %T then
    // если надо, рисуем линии параметров случайного процесса
    mx = mean(x);
    sigmaX = sqrt(corr(x, 0));
    plot2d([1 : n], ones(count, 1) * mx, 5);
    plot2d([1 : n], ones(count, 1) * (mx - sigmaX), 2);
    plot2d([1 : n], ones(count, 1) * (mx + sigmaX), 2);
  end;
endfunction

// строим графики реализаций
scf(1); // эта функция открывает новое графическое окно
plotRealization(x, 120, 1, %T);
plotRealization(imitAR, 120, 15, %F);
// эта функция подписывает цвета линий на графике
legends(['Source         ', ..
         'AR(' + string(bestARparams(1)) + ')          ', ..
         'Average        ', ..
         'St. deviation'], ..
        [1, 15, 5, 2], opt='ll');
// эта функция подписывает координатные оси
xtitle('', 'Index number', 'Random sequence values');

scf(2);
plotRealization(x, 120, 1, %T);
plotRealization(imitMA, 120, 15, %F);
legends(['Source         ', ..
         'MA(' + string(bestMAparams(2)) + ')          ', ..
         'Average        ', ..
         'St. deviation'], ..
        [1, 15, 5, 2], opt='ll');
xtitle('', 'Index number', 'Random sequence values');

scf(3);
plotRealization(x, 120, 1, %T);
plotRealization(imitARMA, 120, 15, %F);
legends(['Source         ', ..
         'ARMA(' + string(bestARMAparams(1)) + ', ' + ..
                   string(bestARMAparams(2)) + ')     ', ..
         'Average        ', ..
         'St. deviation'], ..
        [1, 15, 5, 2], opt='ll');
xtitle('', 'Index number', 'Random sequence values');

printf("Графики реализаций построены!\n");


//----------------------------------------------
printf("\n----------- Задание 6 -----------\n");
//----------------------------------------------

// строим графики нормированных корреляционных функций
scf(4);
plotNormCorrFunc(x, 10, 5);
plotTheorNormCorr(alphasAR, betasAR, 10, 2);
plotNormCorrFunc(imitAR, 10, 15);
plot2d([0 : 10], ones(11, 1) / %e, style = 1);
plot2d([0 : 10], ones(11, 1) * (-1) / %e, style = 1);
legends(['Source       ', ..
         'Imitation    ', ..
         'AR(' + string(bestARparams(1)) + ')        ', ..
         'y = (+/-) 1/e'], ..
        [5, 2, 15, 1], opt='ur');
xtitle('', 'Index number', 'Normalized correlation function');

scf(5);
plotNormCorrFunc(x, 10, 5);
plotTheorNormCorr(alphasMA, betasMA, 10, 2);
plotNormCorrFunc(imitMA, 10, 15);
plot2d([0 : 10], ones(11, 1) / %e, style = 1);
plot2d([0 : 10], ones(11, 1) * (-1) / %e, style = 1);
legends(['Source       ', ..
         'Imitation    ', ..
         'MA(' + string(bestMAparams(2)) + ')        ', ..
         'y = (+/-) 1/e'], ..
        [5, 2, 15, 1], opt='ur');
xtitle('', 'Index number', 'Normalized correlation function');

scf(6);
plotNormCorrFunc(x, 10, 5);
plotTheorNormCorr(alphasARMA, betasARMA, 10, 2);
plotNormCorrFunc(imitARMA, 10, 15);
plot2d([0 : 10], ones(11, 1) / %e, style = 1);
plot2d([0 : 10], ones(11, 1) * (-1) / %e, style = 1);
legends(['Source       ', ..
         'Imitation    ', ..
         'ARMA(' + string(bestARMAparams(1)) + ', ' + ..
                   string(bestARMAparams(2)) + ') ', ..
         'y = (+/-) 1/e'], ..
        [5, 2, 15, 1], opt='ur');
xtitle('', 'Index number', 'Normalized correlation function');

printf("Графики нормированных корреляционных функций построены!\n");

// строим графики спектральных плотностей мощности
scf(7);
plotSpecPowDens(x, 10, 0.01, 5); // красный, исходная выборка
plotARMASpecPowDens(alphasAR, betasAR, 0.01, 15); // зеленый, модель АРСС
plotSpecPowDens(imitAR, 10, 0.01, 2); // синий, сгенерированная выборка
legends(['Source    ', ..
         'AR(' + string(bestARparams(1)) + ')     ', ..
         'Imitation'], ..
        [5, 15, 2], opt='ur');
xtitle('', 'Frequency', 'Normalized power spectral density');

scf(8);
plotSpecPowDens(x, 10, 0.01, 5); // красный, исходная выборка
plotARMASpecPowDens(alphasMA, betasMA, 0.01, 15); // зеленый, модель АРСС
plotSpecPowDens(imitMA, 10, 0.01, 2); // синий, сгенерированная выборка
legends(['Source    ', ..
         'MA(' + string(bestMAparams(2)) + ')     ', ..
         'Imitation'], ..
        [5, 15, 2], opt='ur');
xtitle('', 'Frequency', 'Normalized power spectral density');

scf(9);
plotSpecPowDens(x, 10, 0.01, 5); // красный, исходная выборка
plotARMASpecPowDens(alphasARMA, betasARMA, 0.01, 15); // зеленый, модель АРСС
plotSpecPowDens(imitARMA, 10, 0.01, 2); // синий, сгенерированная выборка
legends(['Source       ', ..
         'ARMA(' + string(bestARMAparams(1)) + ', ' + ..
                   string(bestARMAparams(2)) + ') ', ..
         'Imitation'], ..
        [5, 15, 2], opt='ur');
xtitle('', 'Frequency', 'Normalized power spectral density');

printf("Графики спектральных плотностей мощности построены!\n");

// вычисляем все необходимые значения
// минимумы
minX = min(x);
minAR = min(imitAR);
minMA = min(imitMA);
minARMA = min(imitARMA);
// максимумы
maxX = max(x);
maxAR = max(imitAR);
maxMA = max(imitMA);
maxARMA = max(imitARMA);
// средние значения
mAR = mean(imitAR);
mMA = mean(imitMA);
mARMA = mean(imitARMA);
// дисперсии
dAR = corr(imitAR, 0);
dMA = corr(imitMA, 0);
dARMA = corr(imitARMA, 0);
// сред. квадр. отклонения
sigmaX = sqrt(dx);
sigmaAR = sqrt(dAR);
sigmaMA = sqrt(dMA);
sigmaARMA = sqrt(dARMA);
// нормированные корреляционные функции
r = getNormCorrValues(x, 10);
r_AR = getNormCorrValues(imitAR, 10);
r_MA = getNormCorrValues(imitMA, 10);
r_ARMA = getNormCorrValues(imitARMA, 10);
theor_r_AR = theorNormCorr(alphasAR, betasAR, 10);
theor_r_MA = theorNormCorr(alphasMA, betasMA, 10);
theor_r_ARMA = theorNormCorr(alphasARMA, betasARMA, 10);

// и выводим их на экран
printfVector("\nМинимумы (исходный процесс, АР, СС, АРСС):\n", ..
             [minX, minAR, minMA, minARMA]);
printfVector("\nМаксимумы (исходный процесс, АР, СС, АРСС):\n", ..
             [maxX, maxAR, maxMA, maxARMA]);
printfVector("\nСредние значения (исходный процесс, АР, СС, АРСС):\n", ..
             [mx, mAR, mMA, mARMA]);
printfVector("\nДисперсии (исходный процесс, АР, СС, АРСС):\n", ..
             [dx, dAR, dMA, dARMA]);
printfVector("\nСред. квадр. отклонения " + ..
             "(исходный процесс, АР, СС, АРСС):\n", ..
             [sigmaX, sigmaAR, sigmaMA, sigmaARMA]);
 
printfMatrix("\nЗначения нормированных корреляционных функций:\n" + ..
             "(Исходный процесс, теория АР, выборка АР, " + ..
             "теория СС, выборка СС, теория АРСС, выборка АРСС)", ..
[r, theor_r_AR, r_AR, theor_r_MA, r_MA, theor_r_ARMA, r_ARMA]);

// СКО корреляционных функций
deviations = [deviation(r, r, 10), ..
              deviation(r, theor_r_AR, 10), ..
              deviation(r, r_AR, 10), ..
              deviation(r, theor_r_MA, 10), ..
              deviation(r, r_MA, 10), ..
              deviation(r, theor_r_ARMA, 10), ..
              deviation(r, r_ARMA, 10)
             ];
printfVector("\nСКО корреляционных функций:\n", deviations);

//----------------------------------------------
printf("\n------ Задание 6 выполнено ------\n");
//----------------------------------------------

