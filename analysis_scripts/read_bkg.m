function bkg = read_bkg(filename)

error(nargchk(1,1,nargin,'struct'));

[fn, basepath, suffix, exists] = check_filename(filename, 'ascbkg');

if exists==0
  bkg = struct([]);
  warning(['read_bkg: File ' filename ' does not exist']);
  return
end

fid = fopen(fn, 'r');
bkg.type = 'ascbkg';
try
  temp = fscanf(fid, '%f', 3);

  bkg.nSHOT = temp(1);
  bkg.tSHOT = temp(2);
  bkg.modflg = temp(3);

  fgetl(fid);
  bkg.devnam = fgetl(fid);
  
  temp = fscanf(fid, '%f', 2);

  bkg.FPPkat = temp(1);
  bkg.IpiFPP = temp(2);

  LPFx = fscanf(fid, '%d', 1);
  
  bkg.PFxx = fscanf(fid, '%f', [1 LPFx+1]);
  bkg.RPFx = fscanf(fid, '%f', [1 LPFx+1]);
  bkg.zPFx = fscanf(fid, '%f', [1 LPFx+1]);

  bkg.SSQ = fscanf(fid, '%f', [1 4]);

  LPF = fscanf(fid, '%d', 1);
  bkg.rhoPF = fscanf(fid, '%f', [LPF, 1]);
  bkg.PFL = fscanf(fid, '%f', [LPF, 1]);
  bkg.Ne = fscanf(fid, '%f', [LPF, 1]);
  bkg.Te = fscanf(fid, '%f', [LPF, 1]);
  bkg.Ni = fscanf(fid, '%f', [LPF, 1]);
  bkg.Ti = fscanf(fid, '%f', [LPF, 1]);
  bkg.Vol = fscanf(fid, '%f', [LPF, 1]);
  bkg.Area = fscanf(fid, '%f', [LPF, 1]);
  bkg.Qpl = fscanf(fid, '%f', [LPF, 1]);

  temp = fscanf(fid, '%d', 2);
  linr = temp(1);
  linz = temp(2);

  bkg.tr = fscanf(fid, '%f', [linr 1]);
  bkg.tz = fscanf(fid, '%f', [linz 1]);
  bkg.tbr = fscanf(fid, '%f', [linz linr]);
  bkg.tbz = fscanf(fid, '%f', [linz linr]);
  bkg.tbt = fscanf(fid, '%f', [linz linr]);
  bkg.trhoRZ = fscanf(fid, '%f', [linz linr]);
  bkg.tpsiRZ = fscanf(fid, '%f', [linz linr]);
catch
  warning on
  warning(sprintf('%s\n%s', ...
                  ['Error occurred while reading ' fn ':'], ...
                  ['         ' lasterr]));
end

fclose(fid);