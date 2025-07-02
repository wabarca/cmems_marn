%% Copernicus Marine products
% Por William Abarca, 2025
%% Cleaning
clearvars; close all; clc

%% Area
% --minimum-longitude -110 --maximum-longitude -80 --minimum-latitude 5 --maximum-latitude 20
minLon = -110;
maxLon = -80;
minLat = 5;
maxLat = 20;

%% Dataset Id
% Corrientes, salinidad y SST: cmems_obs-sst_glo_phy_nrt_l4_P1D-m
% Ocean color: cmems_obs-oc_glo_bgc-plankton_nrt_l4-gapfree-multi-4km_P1D
% Sea surface height: cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D
currents = 'cmems_obs-sst_glo_phy_nrt_l4_P1D-m';
oceancolor = 'cmems_obs-oc_glo_bgc-plankton_nrt_l4-gapfree-multi-4km_P1D';
seasurfaceheight = 'cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D';

product = oceancolor;

%% Date 
year = datetime('now','Format','yyyy');
dateini = datetime('now','Format','yyyy-MM-dd');
dateend = dateini + days(7);

switch product
    case currents
        % Output Dir
        ncpath = ['datasets/currents/sst_glo_phy_nrt_l4_P1D/',char(year),'/',char(dateini),'/'];
        copernicus = ['copernicusmarine.exe subset --credentials-file '.copernicusmarine-credentials' --dataset-id '];' ...
            'copernicusmarine subset --credentials-file ".copernicusmarine-credentials" --dataset-id cmems_obs-sst_glo_phy_nrt_l4_P1D-m --variable analysed_sst --start-datetime 2024-01-01T00:00:00 --end-datetime 2024-12-31T00:00:00 --minimum-longitude -92 --maximum-longitude -88 --minimum-latitude 10 --maximum-latitude 14
    case oceancolor
        copernicusmarine subset --credentials-file ".copernicusmarine-credentials" --dataset-id cmems_obs-oc_glo_bgc-plankton_nrt_l4-gapfree-multi-4km_P1D --variable CHL --variable CHL_uncertainty --variable flags --start-datetime 2025-02-18T00:00:00 --end-datetime 2025-02-18T00:00:00 --minimum-longitude -120 --maximum-longitude -60 --minimum-latitude 0 --maximum-latitude 30
    case seasurfaceheight
        copernicusmarine subset --credentials-file ".copernicusmarine-credentials" --dataset-id cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D --variable ugos --variable vgos --variable sla --variable adt --start-datetime 2025-02-19T00:00:00 --end-datetime 2025-02-19T00:00:00 --minimum-longitude -120 --maximum-longitude -60 --minimum-latitude 0 --maximum-latitude 30
    otherwise
        disp('Seleccione un producto valido')
end






copernicus=['copernicusmarine.exe subset --credentials-file '.copernicusmarine-credentials' --dataset-id cmems_mod_glo_phy_anfc_0.083deg_PT1H-m --variable so --variable thetao --variable uo --variable vo --variable zos --start-datetime 2025-02-13T12:00:00 --end-datetime 2025-02-20T12:00:00 --minimum-longitude -110 --maximum-longitude -80 --minimum-latitude 5 --maximum-latitude 20']; 
system(copernicus);
