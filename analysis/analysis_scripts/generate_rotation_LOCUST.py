#generate rotation profile

import pyuda, pathlib, context
import numpy as np
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.rotation import Rotation
import processing.utils

filename_eq = pathlib.Path("LOCUST") / "29034W03" / "g29034"
eq = Equilibrium("", data_format="GEQDSK", filename=filename_eq, GEQDSKFIX=1)


client = pyuda.Client()
shot = 29034
time_slice = 0.36
rot_data = client.get("ACT_SS_VELOCITY", shot)
rot_data_R = rot_data.dims[1].data
rot_time = rot_data.time.data
time_slice_index = np.argmin(np.abs(rot_time - time_slice))
rot = rot_data.data[time_slice_index]

rot_ang = rot / rot_data_R  # omega = 2pi/T = V/R_major

# calculate poloidal flux along rot_data_R at Z=0 (where SS beam is pointed and the CXRS diagnostic looks down)
flux_pol = processing.utils.value_at_RZ(
    R=rot_data_R, Z=np.full(len(rot_data_R), 0.0), quantity=eq["psirz"], grid=eq
)
flux_pol_norm = (flux_pol - eq["simag"]) / (eq["sibry"] - eq["simag"])
r = Rotation("")
r.set(
    rotation_vel=rot,
    rotation_ang=rot_ang,
    R_1D=rot_data_R,
    flux_pol=flux_pol,
    flux_pol_norm=flux_pol_norm,
)

(
    r["flux_pol_norm"],
    r["flux_pol"],
    r["R_1D"],
    r["rotation_vel"],
    r["rotation_ang"],
) = processing.utils.sort_arrays(
    r["flux_pol_norm"], r["flux_pol"], r["R_1D"], r["rotation_vel"], r["rotation_ang"]
)

r.plot()

r.dump_data('LOCUST','profile_wT.dat')