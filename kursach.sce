
stacksize('max');
clear();
funcprot(0);
printf("\n\n--- (�) Alexei Dergunov, 2010 ---\n");

//----------------------------------------------
printf("\n----------- ������� 1 -----------\n");
//----------------------------------------------

x = fscanfMat('C:\data2.txt'); // ��������� ������� �� �����

n = length(x); // ����� �������
if n <> 5000 then
  printf("������� ������� ������� �� �����!\n");
end;

// ������� �������������� �������
function[ans] = corr(x, t)
  ans = 0;
  k = abs(t); // �.�. �������������� ������� ������, R(-t) = R(t)
  n = length(x);
  mx = mean(x);
  // ������ ������� �� �������
  for i = 1 : n - k
    ans = ans + (x(i) - mx) * (x(i + k) - mx);
  end;
  ans = ans / (n - k);
endfunction

// ������� ������������� �������������� �������
function[ans] = normCorr(x, t)
  n = length(x);
  dx = corr(x, 0);
  // ����� �� ���������
  ans = corr(x, t) / dx;
endfunction

// ������� �������� �������������� �������
// R(t), t = [0, tmax]
function[ans] = getCorrValues(x, tmax)
  ans = zeros(tmax + 1, 1);
  for i = 0 : tmax
    ans(i + 1) = corr(x, i);
  end;
endfunction

// ������� �������� ������������� �������������� �������
// r(t), t = [0, tmax]
function[ans] = getNormCorrValues(x, tmax)
  ans = zeros(tmax + 1, 1);
  for i = 0 : tmax
    ans(i + 1) = normCorr(x, i);
  end;
endfunction

// ������� ������ ����������
function[ans] = corrRadius(x, tmax)
  // �������� ��������� ������ ����������
  // � t = tmax � ������� �������
  for t = tmax : -1 : 0
    if abs(normCorr(x, t)) >= 1 / %e
      // ��� ������ ��������� ��������, ������� 1 / e,
      // ������� �� �����
      ans = t + 1;
      break;
    end;
  end;
endfunction

// ������ ������ ������������� �������������� �������, t = [0, tmax]
function[] = plotNormCorrFunc(x, tmax, _style)
  t = [0 : tmax]; // �������� �� ��� �������
  r = getNormCorrValues(x, tmax); // �������� �� ��� �������
  plot2d(t, r, style = _style);
endfunction

mx = mean(x); // ���������� �������
printf("���������� �������   = %.4f\n", mx);

dx = corr(x, 0); // ���������� ���������
printf("���������� ��������� = %.4f\n", dx);

T = corrRadius(x, 250); // ������ ����������
printf("������ ����������    = %d\n", T);


//----------------------------------------------
printf("\n----------- ������� 2 -----------\n");
//----------------------------------------------

// ������� ��������� beta(1)...beta(M) ������ ����
function[ans] = findBetas(x, M, N)
  if M == 0 then
    ans = [];
  else
    R = zeros(2 * M, 1);
    // ����� ����� ������, �� ������ �������������� ��������
    for i = N - M + 1 : N + M
      // ������ ��� �������� �������������� �������
      R(i + M - N) = corr(x, i);
    end;
    A = zeros(M, M); // ������� �������
    b = zeros(M, 1); // ������� ��������� ������
    for i = 1 : M
      b(i) = -R(M + i);
      for j = 1 : M
        A(i, j) = R(M + i - j);
      end;
    end;
    ans = linsolve(A, b); // ������ ������� Ax + b = 0
  end;
endfunction

// ������� ��������� alpha(0)...alpha(N) ������ ����
// �� ������ 2
function[ans] = findAlphas2(x, M, N, betas)
  
  // ��������������� ������� ��� ������ 2
  // ���� ������������� ������ ��(N)
  // ��� ������������������ y,
  // ��������� �������� ������������������ x
  function[ans] = generateMovingAverage(x, betas)
    n = length(x);
    ans = zeros(n, 1);
    M = length(betas);
    // ������ ������� �� �������
    for i = 1 : n
      ans(i) = x(i);
      j = 1;
      while j <= M & i - j > 0 do
        ans(i) = ans(i) - betas(j) * x(i - j);
        j = j + 1;
      end;
    end;
  endfunction
  
  // ���������� ������������������ y
  y = generateMovingAverage(x, betas);
  
  // ������� �������� �� �������������� ������� R(n), n = [0, N]
  R = getCorrValues(y, N);
  
  // �������, ����������� ��� ������� ���������� �������
  // ��� ������ ���� ����� ���� ...
  function[ans] = funcToZero(alpha)
    ans = zeros(N + 1, 1);
    for i = 0 : N
      ans(i + 1) = -R(i + 1); // ... ������� ���� �������� ���
      for j = 0 : N - i
        ans(i + 1) = ans(i + 1) + alpha(j + 1) * alpha(j + i + 1);
      end;
    end;
  endfunction
  
  // ������ ���������� �������
  [ans, values, info] = fsolve([1 : N + 1], funcToZero);
  if info == 4 then // �������� ����������� �������� - ������� �� ����������
    ans(1) = %nan;
  end;
  
endfunction

// �������� �������������� �������
// ��� ������������� ������ 1
// ��� ������� ����� ��� ���������� alpha
function[ans] = mutCorr(alphas, betas)
  N = length(alphas) - 1;
  M = length(betas);
  ans = zeros(N + 1, 1);
  // ������� �� �������
  for i = 0 : N
    ans(i + 1) = alphas(i + 1);
    m = min(i, M);
    for j = 1 : m
      ans(i + 1) = ans(i + 1) + betas(j) * ans(i - j + 1);
    end;
  end;
endfunction

// ������� ��������� alpha(0)...alpha(N) ������ ����
// �� ������ 1
function[ans] = findAlphas1(x, M, N, betas)
  // ������� �������� �������������� ������� R(m), m = [0, max(M, N)]
  R = getCorrValues(x, max(M, N));

  // �������, ����������� ��� ������� ���������� �������
  // ��� ������ ���� ����� ����
  function[ans] = funcToZero(alphas)
    ans = zeros(N + 1, 1);
    for i = 0 : N
      ans(i + 1) = -R(i + 1); // "��������� ����"
      // ���������, ���������� beta
      for j = 1 : M
        ans(i + 1) = ans(i + 1) + betas(j) * R(abs(j - i) + 1);
      end;
      RR = mutCorr(alphas, betas); // �������� �������������� �������
      // ���������, ���������� alpha
      for j = i : N
        ans(i + 1) = ans(i + 1) + alphas(j + 1) * RR(j - i + 1);
      end;
    end;
  endfunction

  // ������ ���������� �������
  [ans, values, info] = fsolve([1 : N + 1], funcToZero);
  if info == 4 then // �������� ����������� �������� - ������� �� ����������
    ans(1) = %nan;
  end;
  
endfunction

// ��������� ������������ ������
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

// �������� ��������������� ����� ������� � �������
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
        printf("������ �����������\n");
      else
        alphas = findAlphas1(x, M, N, betas);
        if isnan(alphas(1)) == %T then
          printf("������ �� ����������\n");
        else
          printfVector("alphas = ", alphas);
          printfVector("betas  = ", betas);
        end;
      end;
      printf("\n");
    end;
  end;
endfunction

// ������� �� ����� ��������� ���� �������
dispAllModels;


//----------------------------------------------
printf("\n----------- ������� 3 -----------\n");
//----------------------------------------------

// ������������ ������������� �������������� �������,
// ������ �� ������� ��������� ��� - ������
function[ans] = theorCorr(alphas, betas, t)
  M = length(betas);
  N = length(alphas) - 1;
  n = max(M + 1, N + 1); // ���������� ���������
  RR = mutCorr(alphas, betas);
  b = zeros(n, 1); // ������� ��������� ������
  for i = 0 : N
    for j = i : N
      b(i + 1) = b(i + 1) + alphas(j + 1) * RR(j - i + 1);
    end;
  end;
  A = zeros(n, n); // ������� �������
  c = zeros(n, 1); // ������������ ��� R(i) � �������� �������
  c(1) = -1;
  c(2 : M + 1) = betas;
  c(M + 2 : n) = 0;
  // ���������� ������� �������
  for i = 1 : n
    p = i; // ������, ���������� ������
    k = i; // ������, ���������� �����
    for j = 1 : n
      if k > 0 then
      // ���� � c(i) i > 0
        A(i, j) = A(i, j) + c(k);
      end;
      if p <= n & p <> k then
      // ���� � c(i) i <= n
      // � �������, ���������� ����� � ������, �� ���������
        A(i, j) = A(i, j) + c(p);
      end;
      p = p + 1;
      k = k - 1;
    end;
  end;
  R = linsolve(A, b);
  if t + 1 < n then
    // ��� ����������� �������� ��� ��� ��������
    ans = R(1 : t + 1);
  else
    // ���� ��������� ��� �������, �� �������� �� ����� �
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
  dx = ans(1); // ���������
  ans = ans / dx; // ����� �� ���
endfunction

// ������ ������ ������������� �������������
// �������������� �������, t = [0, tmax]
function[] = plotTheorNormCorr(alphas, betas, tmax, _style)
  t = [0 : tmax]; // �������� �� ��� �������
  r = theorNormCorr(alphas, betas, tmax); // �������� �� ��� �������
  plot2d(t, r, style = _style);
endfunction

// ������� �������������� ����������
// ������ m + 1 �������� (�.�. [0, m])
// ������������� �������������� ������� r1 � r2
function[ans] = deviation(r1, r2, m)
  ans = 0;
  for i = 1 : m + 1
    ans = ans + (r1(i) - r2(i)) ^ 2;
  end;
endfunction

// ���� ������ ������ ��, �� � ����
function[EpsMatr, bestARparams, bestMAparams, bestARMAparams] = ..
        findBestModels(Mmax, Nmax)

  // ������� ��� ����������
  EpsMatr = zeros(Mmax + 1, Nmax + 1); // ������� ����������
  for M = 0 : Mmax
    for N = 0 : Nmax
      betas = findBetas(x, M, N);
      if isStable(betas) == %F then
        EpsMatr(M + 1, N + 1) = %inf; // ������ �����������
      else
        alphas = findAlphas1(x, M, N, betas);
        if isnan(alphas(1)) == %T then
          EpsMatr(M + 1, N + 1) = %inf; // ������ �� ����������
        else
          // ���������� ������������� � ����������
          // ������������� �������������� �������
          theor_r = theorNormCorr(alphas, betas, 10);
          r = getNormCorrValues(x, 10);
          EpsMatr(M + 1, N + 1) = deviation(r, theor_r, 10);
        end;
      end;
    end;
  end;

  // ������� ������ ������ ��
  bestARparams = zeros(1, 2); // ��������� M, N ������ ������ ��  
  AReps = EpsMatr(1, 1);
  for M = 1 : Mmax
    if EpsMatr(M + 1, 1) < AReps then
      AReps = EpsMatr(M + 1, 1);
      bestARparams(1) = M;
    end;
  end;

  // ������� ������ ������ ��
  bestMAparams = zeros(1, 2); // ��������� M, N ������ ������ ��
  MAeps = EpsMatr(1, 1);
  for N = 1 : Nmax
    if EpsMatr(1, N + 1) < MAeps then
      MAeps = EpsMatr(1, N + 1);
      bestMAparams(2) = N;
    end;
  end;

  // ������� ������ ������ ����
  bestARMAparams = zeros(1, 2); // ��������� M, N ������ ������ ����
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

// ������� ������ ������
[EpsMatr, bestARparams, bestMAparams, bestARMAparams] = ..
        findBestModels(3, 3);

// �������� ��������������� ����� �������
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

printfMatrix("������������� ����������� ������� ����(M, N):", EpsMatr);
printf("\n��������� ��������� ������ �� (M, N)   = (%d, %d)\n", ..
        bestARparams(1), bestARparams(2));
printf("��������� ��������� ������ CC (M, N)   = (%d, %d)\n", ..
        bestMAparams(1), bestMAparams(2));
printf("��������� ��������� ������ ��CC (M, N) = (%d, %d)\n", ..
        bestARMAparams(1), bestARMAparams(2));

// � ������� �� ���������
// ��
betasAR = ..
  findBetas(x, bestARparams(1), bestARparams(2));
alphasAR = ..
  findAlphas1(x, bestARparams(1), bestARparams(2), betasAR);
// ��
betasMA = ..
  findBetas(x, bestMAparams(1), bestMAparams(2));
alphasMA = ..
  findAlphas1(x, bestMAparams(1), bestMAparams(2), betasMA);
// ����
betasARMA = ..
  findBetas(x, bestARMAparams(1), bestARMAparams(2));
alphasARMA = ..
  findAlphas1(x, bestARMAparams(1), bestARMAparams(2), betasARMA);


//----------------------------------------------
printf("\n----------- ������� 4 -----------\n");
//----------------------------------------------

// ������������ ��������� �������� ������ ����
// � ����������� alphas, betas Fi(exp(iw))
function[ans] = ARMASpecPowDens(alphas, betas, w)
  N = length(alphas) - 1;
  M = length(betas);
  a = 0; // ��������� �����
  for k = 0 : N
    a = a + alphas(k + 1) * exp(%i * w * k);
  end;
  b = 1; // ����������� �����
  for k = 1 : M
    b = b - betas(k) * exp(%i * w * k);
  end;
  // ��� ����� �� ���������
  ans = (abs((a / b)) ^ 2) / theorCorr(alphas, betas, 0);
endfunction

// ������ ������������ ��������� �������� �������,
// ��� ���������� ������������ ������ R �������� ��
function[ans] = specPowDens(R, w)
  n = length(R);
  // ������ ������� �� �������
  ans = R(1);
  for k = 1 : n - 1
    ans = ans + 2 * R(k + 1) * cos(w * k);
  end;
  // ����� �� ���������
  ans = ans / R(1);
endfunction

// �������� ������ �������� ������������ ���������
// �������� ������ ���� � �������������� alphas, betas
// Fi(w), w = [0, pi] � ����� step
function[ans] = getARMASpecPowDensValues(alphas, betas, step)
  n = int(%pi / step) + 1;
  ans = zeros(n, 1);
  w = 0;
  for i = 1 : n
    ans(i) = ARMASpecPowDens(alphas, betas, w);
    w = min(w + step, %pi);
  end;
endfunction

// �������� ������ �������� ������������ ����������
// �������� ������� x Fi(w), w = [0, pi] � ����� step,
// ��������� ������� R(t), t = [0, tmax]
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

// ������ ������ ������������ ����������
// �������� ������ ����
function[] = plotARMASpecPowDens(alphas, betas, step, _style)
  w = [0 : step : %pi];
  Fi = getARMASpecPowDensValues(alphas, betas, step);
  plot2d(w, Fi, style = _style);
endfunction

// ������ ������ ������ ������������
// ���������� �������� ������� x
function[] = plotSpecPowDens(x, tmax, step, _style)
  w = [0 : step : %pi];
  Fi = getSpecPowDensValues(x, tmax, step);
  plot2d(w, Fi, style = _style);
endfunction


//----------------------------------------------
printf("\n----------- ������� 5 -----------\n");
//----------------------------------------------

// ���������� count �������� ���������� �������� ����
function[ans] = generateARMA(alphas, betas, count)
  badCount = 1000;
  n = count + badCount;
  ans = zeros(n, 1);
  M = length(betas);
  N = length(alphas) - 1;
  ksi = grand(n, 1, 'nor', 0, 1); // ksi(i) ~ N(0, 1)
  // ������ ������� �� �������
  for i = 1 : n
    ans(i) = 0;
    // ��������� N + 1 ��������� � alpha(j)
    j = 0;
    while j <= N & i - j > 0 do
      ans(i) = ans(i) + alphas(j + 1) * ksi(i - j);
      j = j + 1;
    end;
    // ��������� M ��������� � beta(j)
    j = 1;
    while j <= M & i - j > 0 do
      ans(i) = ans(i) + betas(j) * ans(i - j);
      j = j + 1;
    end;
  end;
  // ����������� ������ badCount = 1000 ��������
  ans = ans(badCount + 1 : n);
  // � ���������� ���. �������� mx
  ans = ans + mx;
endfunction

// ���������� 5000 �������� ������ ������
imitAR = generateARMA(alphasAR, betasAR, 5000);
imitMA = generateARMA(alphasMA, betasMA, 5000);
imitARMA = generateARMA(alphasARMA, betasARMA, 5000);

// ������ ���������� count �������� ���������� ��������
// ���� params = true, �������� ����� ����� �������� �������� � ���
function[] = plotRealization(x, count, _style, params)
  n = min(count, length(x));
  // ������ ���� ����������
  plot2d([1 : n], x(1 : n), style = _style);
  if params == %T then
    // ���� ����, ������ ����� ���������� ���������� ��������
    mx = mean(x);
    sigmaX = sqrt(corr(x, 0));
    plot2d([1 : n], ones(count, 1) * mx, 5);
    plot2d([1 : n], ones(count, 1) * (mx - sigmaX), 2);
    plot2d([1 : n], ones(count, 1) * (mx + sigmaX), 2);
  end;
endfunction

// ������ ������� ����������
scf(1); // ��� ������� ��������� ����� ����������� ����
plotRealization(x, 120, 1, %T);
plotRealization(imitAR, 120, 15, %F);
// ��� ������� ����������� ����� ����� �� �������
legends(['Source         ', ..
         'AR(' + string(bestARparams(1)) + ')          ', ..
         'Average        ', ..
         'St. deviation'], ..
        [1, 15, 5, 2], opt='ll');
// ��� ������� ����������� ������������ ���
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

printf("������� ���������� ���������!\n");


//----------------------------------------------
printf("\n----------- ������� 6 -----------\n");
//----------------------------------------------

// ������ ������� ������������� �������������� �������
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

printf("������� ������������� �������������� ������� ���������!\n");

// ������ ������� ������������ ���������� ��������
scf(7);
plotSpecPowDens(x, 10, 0.01, 5); // �������, �������� �������
plotARMASpecPowDens(alphasAR, betasAR, 0.01, 15); // �������, ������ ����
plotSpecPowDens(imitAR, 10, 0.01, 2); // �����, ��������������� �������
legends(['Source    ', ..
         'AR(' + string(bestARparams(1)) + ')     ', ..
         'Imitation'], ..
        [5, 15, 2], opt='ur');
xtitle('', 'Frequency', 'Normalized power spectral density');

scf(8);
plotSpecPowDens(x, 10, 0.01, 5); // �������, �������� �������
plotARMASpecPowDens(alphasMA, betasMA, 0.01, 15); // �������, ������ ����
plotSpecPowDens(imitMA, 10, 0.01, 2); // �����, ��������������� �������
legends(['Source    ', ..
         'MA(' + string(bestMAparams(2)) + ')     ', ..
         'Imitation'], ..
        [5, 15, 2], opt='ur');
xtitle('', 'Frequency', 'Normalized power spectral density');

scf(9);
plotSpecPowDens(x, 10, 0.01, 5); // �������, �������� �������
plotARMASpecPowDens(alphasARMA, betasARMA, 0.01, 15); // �������, ������ ����
plotSpecPowDens(imitARMA, 10, 0.01, 2); // �����, ��������������� �������
legends(['Source       ', ..
         'ARMA(' + string(bestARMAparams(1)) + ', ' + ..
                   string(bestARMAparams(2)) + ') ', ..
         'Imitation'], ..
        [5, 15, 2], opt='ur');
xtitle('', 'Frequency', 'Normalized power spectral density');

printf("������� ������������ ���������� �������� ���������!\n");

// ��������� ��� ����������� ��������
// ��������
minX = min(x);
minAR = min(imitAR);
minMA = min(imitMA);
minARMA = min(imitARMA);
// ���������
maxX = max(x);
maxAR = max(imitAR);
maxMA = max(imitMA);
maxARMA = max(imitARMA);
// ������� ��������
mAR = mean(imitAR);
mMA = mean(imitMA);
mARMA = mean(imitARMA);
// ���������
dAR = corr(imitAR, 0);
dMA = corr(imitMA, 0);
dARMA = corr(imitARMA, 0);
// ����. �����. ����������
sigmaX = sqrt(dx);
sigmaAR = sqrt(dAR);
sigmaMA = sqrt(dMA);
sigmaARMA = sqrt(dARMA);
// ������������� �������������� �������
r = getNormCorrValues(x, 10);
r_AR = getNormCorrValues(imitAR, 10);
r_MA = getNormCorrValues(imitMA, 10);
r_ARMA = getNormCorrValues(imitARMA, 10);
theor_r_AR = theorNormCorr(alphasAR, betasAR, 10);
theor_r_MA = theorNormCorr(alphasMA, betasMA, 10);
theor_r_ARMA = theorNormCorr(alphasARMA, betasARMA, 10);

// � ������� �� �� �����
printfVector("\n�������� (�������� �������, ��, ��, ����):\n", ..
             [minX, minAR, minMA, minARMA]);
printfVector("\n��������� (�������� �������, ��, ��, ����):\n", ..
             [maxX, maxAR, maxMA, maxARMA]);
printfVector("\n������� �������� (�������� �������, ��, ��, ����):\n", ..
             [mx, mAR, mMA, mARMA]);
printfVector("\n��������� (�������� �������, ��, ��, ����):\n", ..
             [dx, dAR, dMA, dARMA]);
printfVector("\n����. �����. ���������� " + ..
             "(�������� �������, ��, ��, ����):\n", ..
             [sigmaX, sigmaAR, sigmaMA, sigmaARMA]);
 
printfMatrix("\n�������� ������������� �������������� �������:\n" + ..
             "(�������� �������, ������ ��, ������� ��, " + ..
             "������ ��, ������� ��, ������ ����, ������� ����)", ..
[r, theor_r_AR, r_AR, theor_r_MA, r_MA, theor_r_ARMA, r_ARMA]);

// ��� �������������� �������
deviations = [deviation(r, r, 10), ..
              deviation(r, theor_r_AR, 10), ..
              deviation(r, r_AR, 10), ..
              deviation(r, theor_r_MA, 10), ..
              deviation(r, r_MA, 10), ..
              deviation(r, theor_r_ARMA, 10), ..
              deviation(r, r_ARMA, 10)
             ];
printfVector("\n��� �������������� �������:\n", deviations);

//----------------------------------------------
printf("\n------ ������� 6 ��������� ------\n");
//----------------------------------------------

