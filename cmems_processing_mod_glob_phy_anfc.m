% === Configuración general ===
clearvars; close all; clc
fecha_hoy = datetime('today');  % Hoy
fecha_ini = datestr(fecha_hoy, 'yyyy-mm-dd');
fecha_fin = datestr(fecha_hoy + days(7), 'yyyy-mm-dd');
fecha_str = datestr(fecha_hoy, 'yyyymmdd');

dataset_id = 'cmems_mod_glo_phy_anfc_0.083deg_PT1H-m';
variables = {'so', 'thetao', 'uo', 'vo', 'zos'};
lat_min = 5; lat_max = 20;
lon_min = -110; lon_max = -80;

credenciales = '.copernicusmarine-credentials';

% === Carpeta de destino ===
carpeta_destino = fullfile('Datasets\globphy\', fecha_str, '\');
if ~exist(carpeta_destino, 'dir')
    mkdir(carpeta_destino);
end

% === Verificar si el archivo ya existe y su tamaño ===
nombre_destino = ['cmems_glob_phy_' fecha_str '.nc'];
ruta_destino = fullfile(carpeta_destino, nombre_destino);
descargar = true;

if isfile(ruta_destino)
    info_archivo = dir(ruta_destino);
    if info_archivo.bytes > 125e6  % 125 MB en bytes
        fprintf('Archivo ya existe y supera los 125 MB: %s\n', ruta_destino);
        descargar = false;
    end
end

% === Ejecutar descarga si es necesario ===
if descargar
    % Construir el comando
    cmd = ['copernicusmarine.exe subset ' ...
        '--credentials-file "' credenciales '" ' ...
        '--dataset-id ' dataset_id ' '];

    for i = 1:length(variables)
        cmd = [cmd '--variable ' variables{i} ' '];
    end

    cmd = [cmd ...
        '--start-datetime ' fecha_ini 'T00:00:00 ' ...
        '--end-datetime ' fecha_fin 'T18:00:00 ' ...
        '--minimum-longitude ' num2str(lon_min) ' ' ...
        '--maximum-longitude ' num2str(lon_max) ' ' ...
        '--minimum-latitude ' num2str(lat_min) ' ' ...
        '--maximum-latitude ' num2str(lat_max) ' ' ...
        '--netcdf-compression-level 9'];

    disp('Ejecutando descarga de datos desde Copernicus Marine...');
    [status, result] = system(cmd);
    disp(result);

    if status ~= 0
        error('Error al ejecutar copernicusmarine.exe');
    end

    % Buscar y mover archivo descargado
    archivos = dir('cmems_mod_glo_phy_anfc_0.083deg_PT1H-m_multi-vars_*.nc');
    if isempty(archivos)
        error('No se encontró el archivo descargado.');
    end

    nombre_original = archivos(1).name;
    movefile(nombre_original, ruta_destino);
    disp(['Archivo descargado y movido a: ' ruta_destino]);
else
    disp('Omitiendo descarga. Procediendo al procesamiento...');
end

% === Lectura de variables del archivo NetCDF ===
disp('Leyendo variables del archivo (cada 6 horas)...');

% Leer coordenadas
lon = ncread(ruta_destino, 'longitude');
lat = ncread(ruta_destino, 'latitude');

% Malla de coordenadas
[LON, LAT] = meshgrid(lon,lat);

% Leer tiempo completo y convertir a hora local
time_raw = ncread(ruta_destino, 'time');
fecha_base = datetime(1950,1,1);
time_full = fecha_base + hours(time_raw);
idx_cada_6h = 1:6:length(time_raw);  % índices válidos cada 6h
time = time_full(idx_cada_6h) - hours(6);  % hora local
fecha = datestr(time,'dd-mm-yy HH');

% Leer variables completas (todo el time)
thetao = squeeze(ncread(ruta_destino, 'thetao'));
so  = squeeze(ncread(ruta_destino, 'so'));
uo  = squeeze(ncread(ruta_destino, 'uo'));
vo  = squeeze(ncread(ruta_destino, 'vo'));
zos = squeeze(ncread(ruta_destino, 'zos'));

% Conversión a nudos
mss_a_nudos = 3600 / 1852;

% Submuestreo cada 6h en eje time
sss = so(:, :, idx_cada_6h);   % PSU
sst = thetao(:, :, idx_cada_6h);
u_velocity = uo(:, :, idx_cada_6h) * mss_a_nudos;
v_velocity = vo(:, :, idx_cada_6h) * mss_a_nudos;
ssh = zos(:, :, idx_cada_6h);

% Guardar en estructura
ocean_data.sss = sss;
ocean_data.sst = sst;
ocean_data.u_velocity = u_velocity;
ocean_data.v_velocity = v_velocity;
ocean_data.lon = lon;
ocean_data.lat = lat;
ocean_data.ssh = ssh;
ocean_data.time = time;

% Colormap
[sst_cmap, vel_cmap]= hycom_cmap(100);

% === Nodos de monitoreo
% Coordenadas como [lon_index, lat_index]
idx = [240,237,235,245,242,241,249,247,246,242,252,251,259,256,255,266,264,261];
idy = [105,102,98,102,101,97,101,100,96,100,98,94,98,96,92,98,94,89];

nt = length(ocean_data.time);

% === Magnitud y dirección del vector de velocidad
M_v = sqrt(ocean_data.u_velocity.^2 + ocean_data.v_velocity.^2);  % [time x lat x lon]
Or_v = 90 - atan2d(ocean_data.v_velocity, ocean_data.u_velocity);
Or_v(Or_v < 0) = Or_v(Or_v < 0) + 360;

% === Extraer series de tiempo en cada nodo
N_t = zeros(length(idx), nt);  % temperatura
N_s = zeros(length(idx), nt);  % salinidad
N_v = zeros(length(idx), nt);  % velocidad
N_o = zeros(length(idx), nt);  % dirección

for tt = 1:nt
    for ii = 1:length(idx)
        N_t(ii,tt) = round(ocean_data.sst(idx(ii), idy(ii), tt));
        N_s(ii,tt) = round(ocean_data.sss(idx(ii), idy(ii), tt));
        N_v(ii,tt) = round(M_v(idx(ii), idy(ii), tt),2);
        N_o(ii,tt) = round(Or_v(idx(ii), idy(ii), tt));
    end
end

% === Copiar archivo base GLOPHY.xlsx a carpeta destino
nombre_excel_base = 'GLOPHY.xlsx';
nombre_excel_destino = fullfile(carpeta_destino, 'GLOPHY.xlsx');

if ~isfile(nombre_excel_destino)
    copyfile(nombre_excel_base, nombre_excel_destino);
    fprintf('Archivo base GLOPHY.xlsx copiado a: %s\n', nombre_excel_destino);
else
    fprintf('Archivo GLOPHY.xlsx ya existe en la carpeta destino. Se escribirá en él.\n');
end

% === Crear Excel en la copia de destino
warning('off','MATLAB:xlswrite:AddSheet');
fprintf(1, '%s\n', 'Generando series de tiempo en Excel...');

ws = 0;
fecha_strs = cellstr(datestr(ocean_data.time, 'yyyy-mm-dd HH:MM'));

for ii = 1:3:length(idx)
    ws = ws + 1;
    index = 1:nt;

    % Nodo A
    xlswrite(nombre_excel_destino, {['Nodo ',num2str(ii),': '], ...
        num2str(round(ocean_data.lon(idx(ii)),2)), ...
        num2str(round(ocean_data.lat(idy(ii)),2))}, ws, 'A1');
    xlswrite(nombre_excel_destino, {'FECHA - HORA','TEMP ','SAL ','VEL ','DIR '}, ws, 'A2');
    xlswrite(nombre_excel_destino, fecha_strs(index), ws, 'A3');
    xlswrite(nombre_excel_destino, [N_t(ii,index)', N_s(ii,index)', N_v(ii,index)', N_o(ii,index)'], ws, 'B3');

    % Nodo B
    xlswrite(nombre_excel_destino, {['Nodo ',num2str(ii+1),': '], ...
        num2str(round(ocean_data.lon(idx(ii+1)),2)), ...
        num2str(round(ocean_data.lat(idy(ii+1)),2))}, ws, 'G1');
    xlswrite(nombre_excel_destino, {'FECHA - HORA','TEMP ','SAL ','VEL ','DIR '}, ws, 'G2');
    xlswrite(nombre_excel_destino, fecha_strs(index), ws, 'G3');
    xlswrite(nombre_excel_destino, [N_t(ii+1,index)', N_s(ii+1,index)', N_v(ii+1,index)', N_o(ii+1,index)'], ws, 'H3');

    % Nodo C
    xlswrite(nombre_excel_destino, {['Nodo ',num2str(ii+2),': '], ...
        num2str(round(ocean_data.lon(idx(ii+2)),2)), ...
        num2str(round(ocean_data.lat(idy(ii+2)),2))}, ws, 'M1');
    xlswrite(nombre_excel_destino, {'FECHA - HORA','TEMP ','SAL ','VEL ','DIR '}, ws, 'M2');
    xlswrite(nombre_excel_destino, fecha_strs(index), ws, 'M3');
    xlswrite(nombre_excel_destino, [N_t(ii+2,index)', N_s(ii+2,index)', N_v(ii+2,index)', N_o(ii+2,index)'], ws, 'N3');
end

fprintf(1,'%s\n','Excel finalizado');

% === Generación de gráficas
fprintf(1,'%s\n','Generando graficas de GLOPHY');
%Linea de costa
if exist('landmask.mat','file')
    load('landmask.mat');
end
% Limites politicos
if exist('political_boundaries.mat','file')
    load('political_boundaries.mat');
end
% Tomando las lineas dentro del area de interes
M=[];
N=[];
for i=1:1:length(pb)
    if (pb(i,1)>lon_min) & (pb(i,1)<lon_max) & (pb(i,2)>lat_min) & (pb(i,2)<lat_max)
        M(i,1)=pb(i,1);
        N(i,1)=pb(i,2);
    else
        M(i,1)=NaN;
        N(i,1)=NaN;
    end
end

% === Grafica de Temperatura en la superficie del mar
fprintf(1,'%s\n','Generando graficas SST');
parfor kk = 1:length(time)
    % Extraer datos de temperatura superficial del mar
    SST = sst(:,:,kk)';

    % Crear figura
    hfig1 = figure('Visible', 'off'); % evita abrir la ventana durante el proceso paralelo

    % Generar gráfica
    colormap(sst_cmap);
    pcolor(LON, LAT, SST);
    shading flat;
    clim([25 35]);
    hc = colorbar('Location','SouthOutside');
    title({'Temperatura en la superficie del mar (°C)', fecha(kk,:)});
    xlabel('Longitud');
    ylabel('Latitud');
    axis tight;

    % Cargar límites políticos y costa
    hold on;
    for ii = 1:Nb
        if bounds(ii).level == 1 || bounds(ii).level == 3
            patch(bounds(ii).x, bounds(ii).y, [1 1 1]);
        else
            patch(bounds(ii).x, bounds(ii).y, 'b');
        end
    end
    line(M(:,1), N(:,1), 'color', 'k');

    % Guardar imagen
    fname = ['SST-', fecha(kk,1:11)];
    saveas(hfig1, fullfile(carpeta_destino, fname), 'png');
    close(hfig1);
end

% === Grafica de Salinidad en la superficie del mar
fprintf(1,'%s\n','Generando graficas SSS');
parfor kk = 1:length(time)
    % Extraer datos de salinidad superficial del mar
    SSS = sss(:,:,kk)';

    % Crear figura
    hfig2 = figure('Visible', 'off'); % evita abrir la ventana durante el proceso paralelo

    % Generar gráfica
    colormap(sst_cmap);
    pcolor(LON, LAT, SSS);
    shading flat;
    clim([30 40]);
    hc = colorbar('Location','SouthOutside');
    title({'Salinidad en la superficie del mar (psu)', fecha(kk,:)});
    xlabel('Longitud');
    ylabel('Latitud');
    axis tight;

    % Cargar límites políticos y costa
    hold on;
    for ii = 1:Nb
        if bounds(ii).level == 1 || bounds(ii).level == 3
            patch(bounds(ii).x, bounds(ii).y, [1 1 1]);
        else
            patch(bounds(ii).x, bounds(ii).y, 'b');
        end
    end
    line(M(:,1), N(:,1), 'color', 'k');

    % Guardar imagen
    fname = ['SSS-', fecha(kk,1:11)];
    saveas(hfig2, fullfile(carpeta_destino, fname), 'png');
    close(hfig2);
end
% === Grafica del vector de velocidad de corrientes
fprintf(1,'%s\n','Generando graficas Vect_Vel');
parfor kk = 1:length(time)
    % Extraer datos de velocidad de corrientes
    VEL = M_v(:,:,kk)';
    v_v = v_velocity(:,:,kk)';
    u_v = u_velocity(:,:,kk)';

    % Crear figura
    hfig3 = figure('Visible', 'off'); % evita abrir ventana en proceso paralelo

    % Generar gráfica
    colormap(vel_cmap);
    pcolor(LON, LAT, VEL);
    shading flat;
    hold on;
    quiver(LON(1:6:end,1:6:end), LAT(1:6:end,1:6:end), u_v(1:6:end,1:6:end), v_v(1:6:end,1:6:end), 'color','k', 'AutoScaleFactor', 1.8);
    clim([0 2]);
    hc = colorbar('Location','SouthOutside');
    title({'Velocidad de corrientes en la superficie del mar (nudos)', fecha(kk,:)});
    xlabel('Longitud');
    ylabel('Latitud');
    axis tight;
    xlim([min(lon), max(lon)]);
    ylim([min(lat), max(lat)]);

    % Cargar límites políticos y costa
    hold on;
    for ii = 1:Nb
        if bounds(ii).level == 1 || bounds(ii).level == 3
            patch(bounds(ii).x, bounds(ii).y, [1 1 1]);
        else
            patch(bounds(ii).x, bounds(ii).y, 'b');
        end
    end
    line(M(:,1), N(:,1), 'color', 'k');

    % Guardar imagen
    fname = ['Vect_Vel-', fecha(kk,1:11)];
    saveas(hfig3, fullfile(carpeta_destino, fname), 'png');
    close(hfig3);
end
% === Gráfica de la altura de la superficie del mar
fprintf(1,'%s\n','Generando graficas ZOS');
parfor kk = 1:length(time)
    % Extraer datos de altura de la superficie del mar
    ZOS = zos(:,:,kk)';

    % Crear figura
    hfig4 = figure('Visible', 'off'); % evita abrir la ventana durante el proceso paralelo

    % Generar gráfica
    colormap(sst_cmap);
    pcolor(LON, LAT, ZOS);
    shading flat;
    clim([-0.5 0.95]);
    hc = colorbar('Location','SouthOutside');
    title({'Altura de la superficie del mar (m)', fecha(kk,:)});
    xlabel('Longitud');
    ylabel('Latitud');
    axis tight;

    % Cargar límites políticos y costa
    hold on;
    for ii = 1:Nb
        if bounds(ii).level == 1 || bounds(ii).level == 3
            patch(bounds(ii).x, bounds(ii).y, [1 1 1]);
        else
            patch(bounds(ii).x, bounds(ii).y, 'b');
        end
    end
    line(M(:,1), N(:,1), 'color', 'k');

    % Guardar imagen
    fname = ['ZOS-', fecha(kk,1:11)];
    saveas(hfig4, fullfile(carpeta_destino, fname), 'png');
    close(hfig4);
end
% === Creación de animación GIF
png_files = {dir(fullfile(carpeta_destino, '*.png')).name};

%[~, sorted_indices] = sort([png_files.datenum]);
%png_files = {png_files(sorted_indices).name};

file_name2 = ['animacion_' datestr(datetime('today'),'ddmmyy') '.gif'];

% Forever loop
loops=65535;

% Delay time between images
delay = 1;

% Creating the gif
for i=1:length(png_files)
    if strcmpi('gif',png_files{i}(end-2:end))
        [M,c_map]=imread([carpeta_destino,png_files{i}]);
    else
        a=imread([carpeta_destino,png_files{i}]);
        [M,c_map]= rgb2ind(a,256);
    end
    if i==1
        imwrite(M,c_map,[carpeta_destino,file_name2],'gif','LoopCount',loops,'DelayTime',delay)
        % elseif i==length(file_name)
        %     imwrite(M,c_map,[carpeta_destino,file_name2],'gif','WriteMode','append','DelayTime',delay)
    else
        imwrite(M,c_map,[carpeta_destino,file_name2],'gif','WriteMode','append','DelayTime',delay)
    end
end