import context 
import settings,support,processing.utils
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.input_classes.equilibrium import Equilibrium
import matplotlib.pyplot as plt 


bd_filepath=support.dir_input_files / '157418' / 'LOCUST' / 'response' / 'ptcles.dat_wresp'
eq_filepath=support.dir_input_files / '157418' / 'LOCUST' / 'g157418.03000'

eq=Equilibrium(ID='157418 300ms equilibrium',data_format='GEQDSK',filename=eq_filepath,GEQDSKFIX=1)
bd=Beam_Deposition('3D beam depo particle list',filename=bd_filepath,data_format='LOCUST_FO')
bd.set(V_pitch=processing.utils.pitch_calc(bd,equilibria=[eq]))
fig,axs=plt.subplots(3)

bd.plot(style='scatter',weight=False,axes=['X','Y'],fig=fig,ax=axs[0],real_scale=True,LCFS=eq)
bd.plot(style='scatter',weight=False,axes=['R','V_pitch'],fig=fig,ax=axs[1])
bd.plot(style='scatter',weight=False,axes=['E','V_pitch'],fig=fig,ax=axs[2])
plt.show()