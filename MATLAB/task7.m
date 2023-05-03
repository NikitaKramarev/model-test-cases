close all
clear

% Создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink');

% Вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialElements();

% Mомент времени в [s], в который будут производиться расчёты
epoch = 1000; 

% Расчёт положений всех КА в заданный момент времени
constellation.propagateJ2(epoch);

% Расчет параметров низкоорбитальной/высокоорбитальной группы КА
groupNumber = 2;
[satGroupCoordArray, satGroupCount, planeGroupCount, orbitGroupRadius] = findGroupParameters(constellation, groupNumber);

% Построение матрицы связности для выбранной группы КА
reachabilityMatrix = calcReachabilityMatrix(satGroupCoordArray, planeGroupCount, orbitGroupRadius);

% Расчет кратчайшего пути между двумя случайными КА
satInit = ceil(satGroupCount * rand());  
satFinal = ceil(satGroupCount * rand());   
pathSatCoord = calcShortestPath(satGroupCoordArray, reachabilityMatrix, satInit, satFinal);

% 3D-визуализация найденного кратчайшего пути
plotShortestPath(satGroupCoordArray, reachabilityMatrix, pathSatCoord)


function reachabilityMatrix = calcReachabilityMatrix(satCoordArray, planeCount, orbitRadius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ОПИСАНИЕ:
% Функция calcReachabilityMatrix строит (симметричную) матрицу связности 
% для передачи данных между КА из низкоорбитальной/высокоорбитальной группы. 
% Каждый КА имеет связь с двумя соседними по плоскости КА, а также с 
% ближайшим аппаратом из другой плоскости, т.е. всего устанавливается три связи. 
%
% ВХОДНЫЕ ДАННЫЕ:
% satCoordArray --- декартовы координаты всех КА в [m]. Представляет собой 
% массив размером "satTotalCount" на 3, где "satTotalCount" есть полное 
% число КА на одинаковой высоте над Землей.
% planeCount --- число орбитальных плоскостей (скаляр).
% orbitRadius --- радиус орбиты КА из данного эшелона в [m] (скаляр).
%
% ВЫХОДНЫЕ ЗНАЧЕНИЯ:
% reachabilityMatrix --- матрица связности для межспутниковой связи.
% размером "satTotalCount" на "satTotalCount". Состоит из логических нулей
% (если связи между двумя КА нет) и единиц (если связь есть).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Расчет числа КА в каждой плоскости
    satTotalCount = size(satCoordArray,1);                       % полное число КА
    if (mod(satTotalCount,planeCount) ~= 0)
        error(['Ошибка из функции "calcReachabilityMatrix":' ...
               'число КА в каждой плоскости не равно друг другу!']);
    end
    satPerPlane = satTotalCount / planeCount;                    % число КА в каждой плоскости
    satPerPlaneMinus1 = satPerPlane - 1;

    % Инициализация матрицы связности
    reachabilityMatrix = false(satTotalCount);
    for satIdx = 1:satTotalCount
        reachabilityMatrix(satIdx, satIdx) = true;               % каждый элемент связан сам с собой
    end

    % Установка связи между соседними КА в каждой плоскости (пункт а) задания 7)
    for planeIdx = 1:planeCount
        planeIdxMinus1 = planeIdx - 1;
        for satInPlaneIdx = 1:satPerPlane
            satIdx = satInPlaneIdx + planeIdxMinus1 * satPerPlane;
            if (satInPlaneIdx == satPerPlane)
                satIdxPlus1 = satIdx - satPerPlaneMinus1;
            else
                satIdxPlus1 = satIdx + 1;
            end
            reachabilityMatrix(satIdx,satIdxPlus1) = true;
            reachabilityMatrix(satIdxPlus1,satIdx) = true;
        end % Конец цикла по КА в конкретной плоскости
    end % Конец цикла по плоскостям

    % Проверка выполнения пункта а) задания 7
    for satIdx = 1:satTotalCount        
        connectionCount = (sum(reachabilityMatrix(satIdx, :)) - 1);
        if (connectionCount ~= 2)
            error(['Ошибка из функции "calcReachabilityMatrix" : ' ...
                   'в пределах плоскости связь между соседними КА не установлена!']);
        end
    end

    % Создание матрицы расстояний между КА
    lengthArcMax = pi;                                                      % Максимальное значение углового расстояния между КА
    lengthArcArray = zeros(satTotalCount);                                  % Инициализация матрицы расстояний
    for planeIdx = 1:planeCount                                             
        planeIdxMinus1 = planeIdx - 1;
        satFirst = planeIdxMinus1 * satPerPlane + 1;                        
        satLast = planeIdx * satPerPlane;
        for satIdx1 = satFirst:satLast   
            coordSat1 = satCoordArray(satIdx1, :);                          % Координаты КА-1
            for satIdx2 = 1:satTotalCount                                      
                isDiffer = or(satIdx2 < satFirst, satIdx2 > satLast);
                if isDiffer                                                 % если КА принадлежат разным плоскостям
                    coordSat2 = satCoordArray(satIdx2, :);                  % Координаты КА-2 
                    coordDifferenceOverRadius = (coordSat1 - coordSat2) / orbitRadius;
                    lengthArc = findLengthArc(coordDifferenceOverRadius);   % угловая длина дуги между КА-1 и КА-2
                    lengthArcArray(satIdx1, satIdx2) = lengthArc;
                    lengthArcArray(satIdx2, satIdx1) = lengthArcMax;
                else                                                        % если КА принадлежат одной плоскости
                    lengthArcArray(satIdx1, satIdx2) = lengthArcMax;
                    lengthArcArray(satIdx2, satIdx1) = lengthArcMax;
                end
            end % Конец цикла по всем КА
        end % Конец цикла по КА в конкретной плоскости
    end % Конец цикла по плоскостям

    % Корректировка матрицы связности с учетом межплоскостных связей между КА (пункт б) задания 7)
    lengthArcMin = min(lengthArcArray, [], "all");
    while (lengthArcMin < lengthArcMax)
        [stringMin, columnMin] = find(lengthArcArray == lengthArcMin); 
        satIdx1 = stringMin(1);                           
        satIdx2 = columnMin(1);

        reachabilityMatrix(satIdx1, satIdx2) = true; % устанавливаем связь между КА с минимальным расстоянием 
        reachabilityMatrix(satIdx2, satIdx1) = true;

        for satIdx = 1:satTotalCount
            lengthArcArray(satIdx1, satIdx) = lengthArcMax;
            lengthArcArray(satIdx, satIdx1) = lengthArcMax;
            lengthArcArray(satIdx2, satIdx) = lengthArcMax;
            lengthArcArray(satIdx, satIdx2) = lengthArcMax;
        end
        lengthArcMin = min(lengthArcArray, [], "all");
    end

    % Проверка выполнения пункта б) задания 7
    for satIdx = 1:satTotalCount        
        connectionCount = (sum(reachabilityMatrix(satIdx, :)) - 1);
        if (connectionCount ~= 3)
            error(['Ошибка из функции "calcReachabilityMatrix" : ' ...
                   'не удалось установить связь между КА из разных плоскостей!']);
        end
    end

end



function satPathCoord = calcShortestPath(satCoordArray, reachabilityMatrix, satInit, satFinal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ОПИСАНИЕ:
% Функция calcShortestPath находит координаты КА, через которые проходит 
% кратчайший путь между двумя произвольными КА. Поиск кратчайшего пути 
% осуществляется путем перебора ближайших непомеченных КА в соответствие с
% алгоритмом Дейкстры для равновзвешенного неориентированного графа, 
% имитирующего сеть межспутниковых связей. Подробности см., например, в
% E. W. Dijkstra "A note on two problems in connexion with graphs",
% Numerische Mathematik, volume 1, pages 269–271 (1959).
%
% ВХОДНЫЕ ДАННЫЕ:
% satCoordArray --- декартовы координаты всех КА в [m]. Представляет собой 
% массив размером "satTotalCount" на 3, где "satTotalCount" есть полное 
% число КА на одинаковой высоте над Землей.
% reachabilityMatrix --- матрица связности для межспутниковой связи
% размером "satTotalCount" на "satTotalCount". Состоит из логических нулей
% (если связи между двумя КА нет) и единиц (если связь есть).
% satInit, satFinal --- номера начального и конечного КА, между которыми 
% ищется кратчайший путь.
%
% ВЫХОДНЫЕ ЗНАЧЕНИЯ:
% satPathCoord --- декартовы координаты всех КА в [m], через которые
% проходит кратчайший путь. Представляет собой массив размером 
% "satPathCount" на 3, где "satPathCount" есть полное число КА, составляющих
% искомый кратчайший путь.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Проверка входных номеров спутников
    satTotalCount = size(satCoordArray,1); 
    isEqual = (satInit == satFinal) ;
    if isEqual
        error(['Ошибка из функции "calcShortestPath" : ' ...
               'номера выбранных КА совпадают. Задайте номера КА вручную или просто перезапустите программу.']);
    end
    isOutOfRange1 = ((satInit < 1) | (satInit > satTotalCount));
    isOutOfRange2 = ((satFinal < 1) | (satFinal > satTotalCount));
    if or(isOutOfRange1, isOutOfRange2)
        error(['Ошибка из функции "calcShortestPath" : ' ...
               'номера выбранных КА должны принимать значения от 1 до ' int2str(satTotalCount) '!']);
    end

    % Инициализация вектора с номерами-расстояниями до исходного КА
    pathShortest = zeros(1, satTotalCount);
    for satIdx = 1:satTotalCount 
        if (satIdx ~= satInit)  
            pathShortest(satIdx) = satTotalCount; % помечаем только начальный КА
        end 
    end

    % Поиск кратчайшего пути посредством перебора ближайших непомеченных КА
    isFound = false;
    for pathInd = 0:satTotalCount
        pathExist = false;
        for satIdx1 = 1:satTotalCount
            isExtreme = (pathInd == pathShortest(satIdx1));                % крайний на данной итерации КА
            if isExtreme            
                for satIdx2 = 1:satTotalCount
                    isDifferent = (satIdx1 ~= satIdx2);                    % это не один и тот же КА
                    isConnected = reachabilityMatrix(satIdx1, satIdx2);    % есть связь между КА
                    isUntagged = (pathShortest(satIdx2) == satTotalCount); % ближайший из непомеченных ранее КА
                    if (isDifferent && isConnected && isUntagged)
                        pathShortest(satIdx2) = pathInd + 1;               % помечаем новый КА
                        pathExist = true;
                        if (satIdx2 == satFinal)
                            pathIndMax = pathInd + 1;
                            isFound = true;                                % кратчайший путь между КА найден
                            break                                              
                        end
                    end
                end % Конец цикла по крайним еще не помеченным КА
            end 
            if isFound
                break                                                      % путь найден, выходим из цикла
            end
        end % Конец цикла по крайним помеченным КА
        if isFound
            break                                                          % путь найден, выходим из цикла
        end
        if (pathExist == false)                                            % когда перебрали все КА, а путь не найден
            error(['Ошибка из функции "calcShortestPath" : ' ...
                   'Пути между заданными КА не существует!']);
        end
    end 

    % Номера и координаты КА, через которые проходит кратчайший путь
    satPathCount = pathIndMax + 1;
    satPathNumber = zeros(satPathCount, 1);
    satPathNumber(satPathCount) = satFinal;
    satPathCoord = zeros(satPathCount, 3);
    satPathCoord(satPathCount, :) = satCoordArray(satFinal, :);
    for pathInd = pathIndMax:-1:1
        pathIndMinus1 = pathInd - 1;
        pathIndPlus1 = pathInd + 1;
        for satIdx = 1:satTotalCount
            isDifferent = (satIdx ~= satPathNumber(pathIndPlus1));                 % это не один и тот же КА
            isConnected = reachabilityMatrix(satIdx, satPathNumber(pathIndPlus1)); % между КА-1 и КА-2 есть связь
            isTagged = (pathShortest(satIdx) == (pathIndMinus1));                  % этот КА был помечен ранее меньшим индексом 
            if (isConnected && isDifferent && isTagged)
                satPathNumber(pathInd) = satIdx;
                satPathCoord(pathInd, :) = satCoordArray(satIdx, :);
                break;
            end
        end
    end
    if (satPathNumber(1) ~= satInit)
        error(['Ошибка из функции "calcShortestPath" : ' ...
               'не удалось выделить КА, через которые проходит кратчайший путь!']);
    end

end

function plotShortestPath(satCoordArray, reachabilityMatrix, satPathCoord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ОПИСАНИЕ:
% Функция plotShortestPath визуализирует искомый кратчайший путь между 
% двумя КА (красная пунктирная линия), а также множество всех связей
% между КА (синие сплошные линии).
%
% ВХОДНЫЕ ДАННЫЕ:
% satCoordArray --- декартовы координаты всех КА в [m]. Представляет собой 
% массив размером "satTotalCount" на 3, где "satTotalCount" есть полное 
% число КА на одинаковой высоте над Землей.
% reachabilityMatrix --- матрица связности для межспутниковой связи
% размером "satTotalCount" на "satTotalCount". Состоит из логических нулей
% (если связи между двумя КА нет) и единиц (если связь есть).
% satPathCoord --- вектор с номерами спутников (в пределах от 1 до 
% "satTotalCount"), через которые проходит искомый кратчайший путь.
%
% РЕЗУЛЬТАТ:
% 3D-график в декартовой системе координат, где
% синие точки - положения всех КА,
% синие сплошные линии - межспутниковые связи на основе матрицы связности,
% красные "звездочки" - положения КА, через которые проходит кратчайший путь,
% красные пунктирные линии - межспутниковые связи, составляющие путь
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    satTotalCount = size(satCoordArray, 1);                      % полное число спутников

    % Визуализация положений всех КА
    figure;
    hold on;
    for satIdx = 1:satTotalCount
        thisCoord = satCoordArray(satIdx, :);                    % Декартовы координаты конкретного КА [m]
        plot3(thisCoord(1), thisCoord(2), thisCoord(3),'b.');    % Визуализация положения конкретного КА 
    end

    % Визуализация связей между КА
    for satIdx1 = 1:satTotalCount
        for satIdx2 = satIdx1:satTotalCount
            isDifferent = (satIdx1 ~= satIdx2);                  % Это не один и тот же КА
            isConnected = reachabilityMatrix(satIdx1, satIdx2);  % Между КА-1 и КА-2 есть связь
            if (isDifferent && isConnected)                      % Проверка наличия связи между КА-1 и КА-2
                thisCoord1 = satCoordArray(satIdx1, :);
                thisCoord2 = satCoordArray(satIdx2, :);
                plot3([thisCoord1(1) thisCoord2(1)], ...         % Визуализация связи между КА-1 и КА-2
                      [thisCoord1(2) thisCoord2(2)], ...
                      [thisCoord1(3) thisCoord2(3)],'b-');
            end 
        end % Конец цикла по КА-2
    end % Конец цикла по КА-1

    % Визуализация кратчайшего пути
    plot3(satPathCoord(:, 1), satPathCoord(:, 2), ...
          satPathCoord(:, 3),'r*--','LineWidth', 1.5);

    % Выравнивание и подписывание осей 
    axis equal
    xlabel('X, m');
    ylabel('Y, m');
    zlabel('Z, m');

    % Название графика
    title('3D-визуализация кратчайшего пути между выбранными КА')

    hold off

end

function lengthArc = findLengthArc(coordDifferenceOverRadius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ОПИСАНИЕ:
% Функция findLengthArc возвращает длину геодезической линии (ортодромии) 
% между двумя точками. Форма поверхности предполагается сферической.
%
% ВХОДНЫЕ ДАННЫЕ:
% coordDifferenceOverRadius --- разность декартовых координат двух точек,
% нормированная на радиус сферы. Представляет собой вектор с тремя 
% компонентами.
%
% ВЫХОДНЫЕ ЗНАЧЕНИЯ:
% lengthArc --- угловая длина ортодромии (скаляр), принимающая значения
% в интервале [0, pi].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Расчет длины дуги
    chordLength = sum(coordDifferenceOverRadius.^2); 
    chordLengthOver2 = chordLength / 2;
    lengthArc = 2 * asin(chordLengthOver2);   

    % Проверка выходного значения
    isOutOfRange = or((lengthArc < 0), (lengthArc > pi));
    if isOutOfRange
        error(['Ошибка из функции "findLengthArc" : ' ...
               'угловая длина ортодромии не принадлежит [0; pi] !']);
    end

end        

function [satGroupCoordArray, satGroupCount, planeGroupCount, orbitGroupRadius] = findGroupParameters(constellation, groupNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ОПИСАНИЕ:
% Функция findGroupParameters выделяет в отдельный массив координаты 
% низкоорбитальной/высокоорбитальной группы КА. Также находятся некоторые 
% важные параметры группы КА для дальнейшего построения матрицы связности 
% и нахождения кратчайшего пути между двумя случайными КА (см. функции
% calcReachabilityMatrix и calcShortestPath).
% 
% ВХОДНЫЕ ДАННЫЕ:
% constellation --- объект с информацией обо всех КА "Starlink".
% groupNumber --- номер группировки КА принимает значения 1 и 2 для 
% низкоорбитальной (высота 500 км) и высокоорбитальной  (720 км) группы,
% соответсвенно.
%
% ВЫХОДНЫЕ ЗНАЧЕНИЯ:
% satGroupCoordArray --- декартовы координаты выбранной группы КА 
% в [m]. Представляет собой массив размером "satGroupCount" на 3.
% satGroupCount --- полное число КА данной группы (скаляр).
% planeGroupCount --- число орбитальных плоскостей спутников (скаляр).
% orbitGroupRadius --- радиус орбиты КА из данного эшелона в [m] (скаляр).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Проверяем корректность выбранного номера группы
    isNot1 = (groupNumber ~= 1);
    isNot2 = (groupNumber ~= 2);
    if (isNot1 && isNot2)
        error(['Ошибка из функции "findGroupParameters" : ' ...
               'номер группы КА должен принимать значения "1" или "2" !']);
    end

    % Считаем параметры выбранной группы КА, используемые при выполнении задания 7
    satGroup1Count = constellation.groups{1}.totalSatCount;                 % Полное число КА в первой группе
    satGroup2Count = constellation.groups{2}.totalSatCount;                 % Полное число КА во второй группе
    planeGroupCount = constellation.groups{groupNumber}.planeCount;         % Число плоскостей
    orbitGroupAltitude = constellation.groups{groupNumber}.altitude * 1000; % Высота над Землей в [m]
    orbitGroupRadius = constellation.earthRadius + orbitGroupAltitude;      % Радиус орбиты КА в [m]

    % Создаем массив с координатами КА из выбранной группы
    if (groupNumber == 1)
        satGroupCount = satGroup1Count;
        satGroupCoordArray = constellation.state.eci(1: satGroup1Count, :); 
    else
        satGroupCount = satGroup2Count;
        satGroupCoordArray = constellation.state.eci(satGroup1Count + 1: constellation.totalSatCount, :); 
    end

end

