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
        if info_archivo.bytes > 10e6  % mínimo 10 MB
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

% ===== Lectura de variables 
