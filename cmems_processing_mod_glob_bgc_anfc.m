% === Configuración general ===
clearvars; close all; clc

% Fechas y dominio
fecha_hoy = datetime('today');  % Hoy
fecha_ini = datestr(fecha_hoy - days(2), 'yyyy-mm-dd');
fecha_fin = datestr(fecha_hoy + days(4), 'yyyy-mm-dd');
fecha_str = datestr(fecha_hoy, 'yyyymmdd');

lat_min = 0; lat_max = 30;
lon_min = -120; lon_max = -60;

credenciales = '.copernicusmarine-credentials';

% === Datasets y variables a descargar ===
descargas = {
    'cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m', 'nppv', 'bio';
    'cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m', 'ph',   'car';
    'cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m', 'spco2','co2';
    };

% === Carpeta de destino base ===
carpeta_destino = fullfile('Datasets\globbgc\', fecha_str, '\');
if ~exist(carpeta_destino, 'dir')
    mkdir(carpeta_destino);
end

% === Iterar sobre cada dataset ===
for i = 1:size(descargas,1)
    dataset_id = descargas{i,1};
    variable   = descargas{i,2};
    tipo       = descargas{i,3};  % bio, car, co2

    nombre_destino = ['cmems_glob_bgc_' tipo '_' fecha_str '.nc'];
    ruta_destino = fullfile(carpeta_destino, nombre_destino);

    descargar = true;
    if isfile(ruta_destino)
        info_archivo = dir(ruta_destino);
        if info_archivo.bytes > 10e4  % mínimo 10 MB
            fprintf('Archivo ya existe y es válido: %s\n', ruta_destino);
            descargar = false;
        end
    end

    if descargar
        % Construir comando
        cmd = ['copernicusmarine.exe subset ' ...
            '--credentials-file "' credenciales '" ' ...
            '--dataset-id ' dataset_id ' ' ...
            '--variable ' variable ' ' ...
            '--start-datetime ' fecha_ini 'T00:00:00 ' ...
            '--end-datetime ' fecha_fin 'T18:00:00 ' ...
            '--minimum-longitude ' num2str(lon_min) ' ' ...
            '--maximum-longitude ' num2str(lon_max) ' ' ...
            '--minimum-latitude ' num2str(lat_min) ' ' ...
            '--maximum-latitude ' num2str(lat_max) ' ' ...
            '--minimum-depth 0.4940253794193268 ' ...
            '--maximum-depth 0.4940253794193268 ' ...
            '--netcdf-compression-level 9'];

        fprintf('\nDescargando %s (%s)...\n', variable, dataset_id);
        [status, result] = system(cmd);
        disp(result)

        if status ~= 0
            warning('Fallo la descarga de %s\n', variable);
            continue;
        end

        % Buscar archivo recién descargado
        patron = sprintf('%s_%s_*.nc', dataset_id, variable);
        archivos = dir(patron);
        if isempty(archivos)
            warning('No se encontró el archivo descargado para %s\n', variable);
            continue;
        end

        % Mover y renombrar
        nombre_original = archivos(1).name;
        movefile(nombre_original, ruta_destino);
        fprintf('Archivo guardado como: %s\n', ruta_destino);
    else
        fprintf('Omitiendo descarga de %s. Archivo ya existe.\n', variable);
    end
end

% Lectura de variables
fprintf('\nLeyendo variables descargadas...\n');

% Coordenadas comunes
archivo_base = fullfile(carpeta_destino, ['cmems_glob_bgc_bio_' fecha_str '.nc']);

lon  = ncread(archivo_base, 'longitude');
lat  = ncread(archivo_base, 'latitude');

ocean_data.lon = lon;
ocean_data.lat = lat;


% Malla de coordenadas
[LON, LAT] = meshgrid(ocean_data.lon,ocean_data.lat);

% Leer tiempo completo y convertir a hora local
time_raw = ncread(archivo_base, 'time');
fecha_base = datetime(1950,1,1);
time_full = fecha_base + hours(time_raw) - hours(6) ;
fecha = datestr(time_full,'dd-mm-yy HH');
ocean_data.time = time_full;

% Leyendo variables de cada archivo
variables = {'nppv', 'ph', 'spco2'};
for i = 1:length(variables)
    var = variables{i};
    tipo = '';  % para encontrar el archivo correcto según variable

    switch var
        case 'nppv'
            tipo = 'bio';
        case 'ph'
            tipo = 'car';
        case 'spco2'
            tipo = 'co2';
    end

    archivo = fullfile(carpeta_destino, ['cmems_glob_bgc_' tipo '_' fecha_str '.nc']);
    if isfile(archivo)
        fprintf('Leyendo %s desde %s\n', var, archivo);
        datos = ncread(archivo, var);  % dims: lon x lat x [depth?] x time

        if ndims(datos) == 4  % eliminar dimensión depth (es 1)
            datos = squeeze(datos);
        end

        ocean_data.(var) = datos;
    else
        warning('Archivo no encontrado: %s', archivo);
    end
end

spco2 = ocean_data.spco2;
ph = ocean_data.ph;
nppv = ocean_data.nppv;

% Nodos de monitoreo
% Coordenadas como [lon_index, lat_index]
idx = [121,120,119,122,121,121,124,123,123,121,125,124,127,126,126,129,129,128];
idy = [55,55,53,54,54,53,54,54,53,54,53,52,53,53,51,53,52,50];

% Extrayendo series de tiempo en cada nodo
nt = length(ocean_data.time);

N_nppv = zeros(length(idx),nt);
N_ph = zeros(length(idx),nt);
N_co2 = zeros(length(idx),nt);

for tt = 1:nt
    for ii = 1:length(idx)
        N_nppv(ii,tt) = round(ocean_data.nppv(idx(ii),idy(ii),tt));
        N_co2(ii,tt) = round(ocean_data.spco2(idx(ii),idy(ii),tt));
        N_ph(ii,tt) = round(ocean_data.ph(idx(ii),idy(ii),tt));
    end
end

% === Copiar archivo base GLOBGC.xlsx a carpeta destino
nombre_excel_base = 'GLOBGC.xlsx';
nombre_excel_destino = fullfile(carpeta_destino, 'GLOBGC.xlsx');

if ~isfile(nombre_excel_destino)
    copyfile(nombre_excel_base, nombre_excel_destino);
    fprintf('Archivo base GLOBGC.xlsx copiado a: %s\n', nombre_excel_destino);
else
    fprintf('Archivo GLOBGC.xlsx ya existe en la carpeta destino. Se escribirá en él.\n');
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
    xlswrite(nombre_excel_destino, {'FECHA - HORA','PH ','SPCO2 ','NPPV ',' '}, ws, 'A2');
    xlswrite(nombre_excel_destino, fecha_strs(index), ws, 'A3');
    xlswrite(nombre_excel_destino, [N_ph(ii,index)', N_co2(ii,index)', N_nppv(ii,index)'], ws, 'B3');

    % Nodo B
    xlswrite(nombre_excel_destino, {['Nodo ',num2str(ii+1),': '], ...
        num2str(round(ocean_data.lon(idx(ii+1)),2)), ...
        num2str(round(ocean_data.lat(idy(ii+1)),2))}, ws, 'G1');
    xlswrite(nombre_excel_destino, {'FECHA - HORA','PH ','SPCO2 ','NPPV ',' '}, ws, 'G2');
    xlswrite(nombre_excel_destino, fecha_strs(index), ws, 'G3');
    xlswrite(nombre_excel_destino, [N_ph(ii+1,index)', N_co2(ii+1,index)', N_nppv(ii+1,index)'], ws, 'H3');

    % Nodo C
    xlswrite(nombre_excel_destino, {['Nodo ',num2str(ii+2),': '], ...
        num2str(round(ocean_data.lon(idx(ii+2)),2)), ...
        num2str(round(ocean_data.lat(idy(ii+2)),2))}, ws, 'M1');
    xlswrite(nombre_excel_destino, {'FECHA - HORA','PH ','SPCO2 ','NPPV ',' '}, ws, 'M2');
    xlswrite(nombre_excel_destino, fecha_strs(index), ws, 'M3');
    xlswrite(nombre_excel_destino, [N_ph(ii+2,index)', N_co2(ii+2,index)', N_nppv(ii+2,index)'], ws, 'N3');
end

fprintf(1,'%s\n','Excel finalizado');

% === Generación de gráficas
fprintf(1,'%s\n','Generando graficas de GLOBGC');
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

load('boundaries.mat')

% Grafica de Presión parcial de CO2 en superficie del mar
fprintf(1,'%s\n','Generando graficas SPCO2');
parfor kk = 1:length(time_full)
    % Extraer datos de Presión parcial de CO2 en superficie
    SPCO2 = spco2(:,:,kk)';

    % Crear figura
    hfig1 = figure('Visible', 'off'); % evita abrir la ventana durante el proceso paralelo

    % Generar gráfica
    colormap(turbo);
    pcolor(LON, LAT, SPCO2);
    shading flat;
    clim([30 60]);
    hc = colorbar('Location','SouthOutside');
    title({'Presión parcial de CO2 en superficie (Pa)', fecha(kk,:)});
    xlabel('Longitud');
    ylabel('Latitud');
    axis tight;

    % Cargar límites políticos y costa
    hold on;
    for i = 1:length(boundaries)
        b = boundaries{i};
        plot(b.Lon, b.Lat, 'k', 'LineWidth', 0.3); hold on;
    end

    % Guardar imagen
    fname = ['SPCO2-', fecha(kk,1:11)];
    exportgraphics(hfig1, fullfile(carpeta_destino, fname + ".png"), 'Resolution', 300);
    close(hfig1);
end

% Grafica de ph superficie del mar
fprintf(1,'%s\n','Generando graficas PH');
parfor kk = 1:length(time_full)
    % Extraer datos de ph en superficie del mar
    PH = ph(:,:,kk)';

    % Crear figura
    hfig1 = figure('Visible', 'off'); % evita abrir la ventana durante el proceso paralelo

    % Generar gráfica
    colormap(hsv);
    pcolor(LON, LAT, PH);
    shading flat;
    clim([7 9]);
    hc = colorbar('Location','SouthOutside');
    title({'{\itph} en la superficie del mar.', fecha(kk,:)});
    xlabel('Longitud');
    ylabel('Latitud');
    axis tight;

    % Cargar límites políticos y costa
    hold on;
    for i = 1:length(boundaries)
        b = boundaries{i};
        plot(b.Lon, b.Lat, 'k', 'LineWidth', 0.3); hold on;
    end
    % Guardar imagen
    fname = ['PH-', fecha(kk,1:11)];
    exportgraphics(hfig1, fullfile(carpeta_destino, fname + ".png"), 'Resolution', 300);
    close(hfig1);
end

% Grafica de Produccion Primaria Total
fprintf(1,'%s\n','Generando graficas NPPV');
parfor kk = 1:length(time_full)
    % Extraer datos de nppv en superficie del mar
    NPPV = nppv(:,:,kk)';

    % Crear figura
    hfig1 = figure('Visible', 'off'); % evita abrir la ventana durante el proceso paralelo

    % Generar gráfica
    colormap(turbo);
    pcolor(LON, LAT, NPPV);
    shading flat;
    clim([0 100]);
    hc = colorbar('Location','SouthOutside');
    title({'Producción Primaria Total (mg/m³/día)', fecha(kk,:)});
    xlabel('Longitud');
    ylabel('Latitud');
    axis tight;

    % Cargar límites políticos y costa
    hold on;
    for i = 1:length(boundaries)
        b = boundaries{i};
        plot(b.Lon, b.Lat, 'k', 'LineWidth', 0.3); hold on;
    end


    % Guardar imagen
    fname = ['NPPV-', fecha(kk,1:11)];
    exportgraphics(hfig1, fullfile(carpeta_destino, fname + ".png"), 'Resolution', 300);
    close(hfig1);
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