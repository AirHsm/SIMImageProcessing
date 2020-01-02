function config = ReadConfig(data_folder_path)
  addpath('./others/jsonlab-1.5')
  config = config_ini(loadjson(fullfile(data_folder_path,'config.json')));
  config.Data.folderpath = data_folder_path;

  datapath = fullfile(data_folder_path, config.Data.Path);
  config.Reconstruction.DataPath = datapath;
  if ~exist(datapath,'dir')
    mkdir(datapath);
  end

  resultpath = fullfile(data_folder_path,config.Reconstruction.ResultPath);
  config.Reconstruction.ResultPath = resultpath;
  if ~exist(resultpath, 'dir')
    mkdir(resultpath);
  end

  config.Data.folderpath = data_folder_path;
