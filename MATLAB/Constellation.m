classdef Constellation < handle

    properties
        totalSatCount = 0;
        groups = {};
        state;

        % константы
        earthRadius = 6378135;           % Экваториальный радиус Земли [m]
        earthGM = 3.986004415e+14;       % Гравитационный параметр Земли [m3/s2]
        earthJ2 = 1.082626e-3;           % Вторая зональная гармоника геопотенциала
    end

    methods

        function this = Constellation(varargin)
            if isempty(varargin)
                return
            end
            this.loadFromConfigFile(varargin{1});
        end

        function loadFromConfigFile(this, code)
            fileName = 'constellationsTest.json';
            str = fileread(fileName);
            data = jsondecode(str);
            dataThis = [];

            for i = 1:length(data)
                if strcmpi(data(i).name, code)
                    dataThis = data(i);
                    break
                end
            end

            if isempty(dataThis)
                disp('Группировка не найдена в файле');
                return
            end

            for i = 1:size(dataThis.Walkers, 1)
                group.inclination = deg2rad(dataThis.Walkers(i, 1));        % наклонение орбитальной плоскости
                group.satsPerPlane = dataThis.Walkers(i, 2);				% число КА в каждой орбитальной плоскости группы
                group.planeCount = dataThis.Walkers(i, 3);					% число орбитальных плоскостей в группе
                group.f = dataThis.Walkers(i, 4);							% фазовый сдвиг по аргументу широты между КА в соседних плоскостях
                group.altitude = dataThis.Walkers(i, 5);					% высота орбиты
                group.maxRaan = deg2rad(dataThis.Walkers(i, 6));            % максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей)
                group.startRaan = deg2rad(dataThis.Walkers(i, 7));			% прямое восхождение восходящего узла для первой плоскости
                group.totalSatCount = group.satsPerPlane * group.planeCount;

                this.groups{length(this.groups) + 1} = group;                
                this.totalSatCount = this.totalSatCount + group.totalSatCount;
            end
        end

        function getInitialElements(this)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ОПИСАНИЕ:
        % Функция getInitialElements рассчитывает элементы орбит всех КА 
        % в начальный момент времени. 
        % 
        % ВХОДНЫЕ ДАННЫЕ:
        % Общие параметры каждой группы КА:
        % inclination --- наклонение орбитальной плоскости,
        % satsPerPlane --- число КА в каждой орбитальной плоскости группы,
        % planeCount --- число орбитальных плоскостей в группе,
        % и другие (см. функцию loadFromConfigFile).
        %
        % ВЫХОДНЫЕ ЗНАЧЕНИЯ:
        % elements --- элементы орбиты группировки КА в ИСО 
        % в начальный момент времени. Представляет собой массив размером
        % [totalSatCount x 6], где totalSatCount есть полное число КА. 
        % Значения столбцов: [радиус орбиты в [m], эксцентриситет орбиты 
        % (равен нулю по умолчанию), (равен нулю по умолчанию), прямое 
        % восхождение восходящего узла орбиты КА, наклон орбиты (в радианах),
        % аргумент широты].
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Инициализация искомого массива с начальными элементами орбиты
            this.state.elements = zeros(this.totalSatCount, 6);
            
            % Перебираем группы КА, расположенные на разных высотах
            shift = 1;
            for group = this.groups
                satInGroupCount = group{1}.totalSatCount;
                satPerPlaneInGroup = group{1}.satsPerPlane;
                ending = shift + satInGroupCount - 1;

                % Создаем массив c прямыми восхождениями восходящих узлов орбит
                raanArray = linspace(group{1}.startRaan, group{1}.startRaan + group{1}.maxRaan, group{1}.planeCount + 1);
                raanArray = mod(raanArray(1:end-1), 2 * pi);

                % Инициализация массива с элементами данной группы
                orbitElementsArray = zeros(satInGroupCount, 6);

                % Расчет радиуса орбиты и аргумента широты
                satInGroupNumber = 1;
                raanIdx = 0;
                for thisRaan = raanArray
                    for satInPlaneNumber = 0:satPerPlaneInGroup-1
                        sma = this.earthRadius + group{1}.altitude * 1000;
                        aol = 2 * pi / satPerPlaneInGroup * satInPlaneNumber + 2 * pi / satInGroupCount * group{1}.f * raanIdx;

                        orbitElementsArray(satInGroupNumber, :) = [sma, 0, 0, thisRaan, group{1}.inclination, aol];
                        satInGroupNumber = satInGroupNumber + 1;
                    end % Конец цикла по КА внутри данной плоскости
                    raanIdx = raanIdx + 1;
                end % Конец цикла по орбитальным плоскостям

                % Заполнение искомого массива элементами данной группы
                this.state.elements(shift:ending,:) = orbitElementsArray;
                shift = ending + 1;
            end % Конец цикла по группам КА
        end        

        function propagateJ2(this, epochs)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ОПИСАНИЕ:
        % Функция propagateJ2 рассчитывает положения КА в заданные моменты 
        % времени (начальный момент времени равен нулю). Произведен 
        % учет несферичности геопотенциала с точностью до второй зональной
        % гармоники J2.
        % 
        % ВХОДНЫЕ ДАННЫЕ:
        % epochs --- моменты времени в [s], в которые рассчитываются
        % положения КА. Представляет собой вектор длиной length(epochs).
        %
        % ВЫХОДНЫЕ ЗНАЧЕНИЯ:
        % state.eci --- декартовы координаты группировки КА в ИСО 
        % в заданные моменты времени. Представляет собой массив размером
        % [totalSatCount x 3 x length(epochs)], где totalSatCount есть
        % полное число КА.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Инициализация массива искомых (декартовых) координат КА
            this.state.eci = zeros(this.totalSatCount, 3, length(epochs));

            % Некоторы элементы невозмущенной орбиты на начальную эпоху
            % (см. функцию getInitialElements)
            sma          = this.state.elements(:, 1);
            inclination  = this.state.elements(:, 5);            
            raan0        = this.state.elements(:, 4);
            aol0         = this.state.elements(:, 6);                          % аргумент широты на начальную эпоху
            meanMotion   = sqrt(this.earthGM ./ sma.^3);                       % среднее движение КA в [1/s]
            correctionJ2 = 1.5 * this.earthJ2 .* (this.earthRadius ./ sma).^2; % поправка за вторую зональную гармонику J2

            % Учет второй зональной гармоники
            raanPrecessionRate = - meanMotion .* correctionJ2 .* cos(inclination);
            draconicOmega      = meanMotion .* (1 - correctionJ2) .* (1 - 4 .* cos(inclination).^2);

            for epochIdx = 1:length(epochs)

                % Расчет аргумента широты и прямого восхождения на заданную
                % эпоху и учет поправки за несферичность геопотенциала
                aol = aol0 + epochs(epochIdx) * draconicOmega;
                raanOmega = raan0 + epochs(epochIdx) * raanPrecessionRate;

                % Переход от орбитальных к относительным координатам
                this.state.eci(:, :, epochIdx)  = [sma .* (cos(aol) .* cos(raanOmega) - sin(aol) .* cos(inclination) .* sin(raanOmega)), ...
                                                   sma .* (cos(aol) .* sin(raanOmega) + sin(aol) .* cos(inclination) .* cos(raanOmega)), ...
                                                   sma .* (sin(aol) .* sin(inclination))];
            end
        end
    end
end
