function [fnout, basepath, suffixout, exists] = ...
    check_filename(fn, suffix)

error(nargchk(1,2,nargin));


switch nargin
 case 1
  if ~isstr(fn)
    error('check_filename: First argument must be a string');
  end
  if isempty(fn)
    error('check_filename: First argument must be non-empty');
  end
  n = nargin;
 case 2
  if isstr(fn)+isstr(suffix)<2
    error('check_filename: Both arguments must be strings');
  end
  if isempty(fn)
    fn = 'this';
  end
  n = nargin - isempty(suffix);
end


switch n
 case 1
  if ~isempty(strfind(fn, '.'))
    [dummy1, dummy2] = strtok(fliplr(fn), '.');
    basepath = fliplr(dummy2(2:end));
    suffixout = check_suffix(fliplr(dummy1));
  else
    basepath = 'this';
    suffixout = check_suffix(fn);
  end
  if exist([basepath '.' suffixout], 'file')==2
    fnout = [basepath '.' suffixout];
    exists = true;
    return
  end
 case 2
  suffixout = check_suffix(fliplr(strtok(fliplr(suffix), '.')));
  if exist([fn '.' suffixout], 'file')==2
    fnout = [fn '.' suffixout];
    basepath = fn;
    exists = true;
    return
  end

  if ~isempty(strfind(fn, '.'))
    [dummy1, dummy2] = strtok(fliplr(fn), '.');
    basepath = fliplr(dummy2(2:end));
    if exist([basepath '.' suffixout], 'file')==2
      fnout = [basepath '.' suffixout];
      exists = true;
      return
    end
  end
end

basepath = '';
suffixout = '';
fnout = '';
exists = false;

return

function suffix = check_suffix(s)

suffarr = strvcat('ascbg2','ascbkg','magn_bkg','magn_header', ...
                  'ascini','inipar','ascinp', 'ascout', ...
                  'ascwall', 'cmd', 'console', 'console_p', 'cx2d_', ...
                  'cxphi', 'cxpit', 'cxpitch', 'cxssl', 'cxtot', ...
                  'cxtot_radial', 'cxtot_radial_slots', 'erad', 'ersc', ...
                  'logfile_p', 'lossflux', 'orbit_loss', 'neutral_loss', ...
                  'ripple_bkg', 'ripple_input', 'sol', 'summary', ...
                  'wall', 'rhodist', 'tdist2d', 'edist', 'lhedist', ...
                  'icedist', 'jdist', 'jrodist', 'jrdist', 'jrtdist', ...
                  'vpodist', 'vdist', 'enedist', 'vpadist', 'vpedist', ...
                  'vdist4d', 'epdist2d', 'pitchdist', 'emdist', ...
                  'ss_vdist', 'ss_vpadist', 'ss_vpedist', 'ss_4d-vdist', ...
                  'fafner', 'birth','testfus', 'rzvdist4d', 'rzpedist4d', ...
                  'plasma_1d', 'plasma_2d', 'endstate', 'particles', ...
                  'orbits', 'wall_2d', 'wall_3d', 'wall_hits', ...
                  'ascwall-3d','magnTest','neutral_1d','neutral_3d',...
                  'andiff','poincarePol','poincareTor',...
                  'debug_orbits','fild','debug_nbilog','q_profile','seed','CDF'); 

I = strmatch(s, suffarr);
if length(I) == 1
  suffix = strtrim(suffarr(I,:));
elseif length(I) == 2
  if s == 'vdist'
    suffix = s;
  end
else
  suffix = '';
end

return
