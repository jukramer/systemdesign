
import math as m

D_fus = 0
L1 = 0
L2 = 0
L3 = 0
S_wexp = 0
S_hexp = 0
S_vexp = 0




Swet_wing = 1.07 * 2 * S_wexp
Swet_fus = ((m.pi * D_fus)/4)(1/(3*L1^2)((4*L1^2+D_fus^2/4)^1.5-(D_fus^3/8))-D_fus+4*L2+2*m.sqrt(L3^2+D_fus^2/4))
Swet_hor = 1.05 * 2 * S_hexp
Swet_vert = 1.05 * 2 * S_vexp



print()


S_wet = Swet_wing + Swet_fus + Swet_hor + Swet_vert
print(S_wet)


