thickness = 0
Sexp = 0 
Atop = 0
Aside = 0
K = 3.4

Swet_wing = Sexp* (1.977+0.53*(thickness))
Swet_fus = K*(Atop +Aside)/2

S_wet = Swet_wing + Swet_fus 