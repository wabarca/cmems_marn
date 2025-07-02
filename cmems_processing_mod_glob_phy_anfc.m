% === Configuración general ===
fecha_hoy = datetime('today');  % Ayer
fecha_ini = datestr(fecha_hoy, 'yyyy-mm-dd');
fecha_fin = datestr(fecha_hoy + days(7), 'yyyy-mm-dd');
fecha_str = datestr(fecha_hoy, 'yyyymmdd');

dataset_id = 'cmems_mod_glo_phy_anfc_0.083deg_PT1H-m';
variables = {'so', 'thetao', 'uo', 'vo', 'zos'};
lat_min = 5; lat_max = 20;
lon_min = -110; lon_max = -80;

credenciales = '.copernicusmarine-credentials';  % ajustar si está en otro directorio

% === Carpeta de destino ===
carpeta_destino = fullfile('Datasets', fecha_str);
if ~exist(carpeta_destino, 'dir')
    mkdir(carpeta_destino);
end

% === Construcción del comando ===
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
    '--netcdf-compression-level ' num2str(9)];

% === Ejecutar descarga ===
disp('Ejecutando descarga de datos desde Copernicus Marine...');
[status, result] = system(cmd);
disp(result);

if status ~= 0
    error('Error al ejecutar copernicusmarine.exe');
end

% === Buscar y mover archivo descargado ===
archivos = dir('cmems_mod_glo_phy_anfc_0.083deg_PT1H-m_multi-vars_*.nc');
if isempty(archivos)
    error('No se encontró el archivo descargado.');
end

nombre_original = archivos(1).name;
nombre_destino = ['cmems_glo_phy_' fecha_str '.nc'];
ruta_destino = fullfile(carpeta_destino, nombre_destino);

movefile(nombre_original, ruta_destino);
disp(['Archivo descargado y movido a: ' ruta_destino]);

% === Listo para procesamiento posterior ===
