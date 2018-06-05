from pygoda import camgoda

# del(Rp) = (Re_bar*del_E + E_bar*del_Re + Rc_bar*del_C + C_bar*del_Rc - Rp_bar*del_P) / P_bar

cfile = 'fc5.2deg.wiso.piControl_kn028/AMWG.climoFiles/fc5.2deg.wiso.piControl_kn028_JAS_climo.nc'
tfile = 'fc5.2deg.wiso.mh6ka_kn032/AMWG.climoFiles/fc5.2deg.wiso.mh6ka_kn032_JAS_climo.nc'

cnc = camgoda(cfile)
tnc = camgoda(tfile)

box = (0, 35, 330, 35)


# Done with "delta values" 
tPRECT_H218O = tnc.variable("PRECT_H218O", box)
tPRECT_H2O  = tnc.variable("PRECT_H2O", box)
tRp = (tPRECT_H218O / tPRECT_H2O - 1) * 1000

tQFLX_H218O = tnc.variable("QFLX_H218O", box)
tQFLX_H2O = tnc.variable("QFLX_H2O", box)
tRe = (tQFLX_H218O / tQFLX_H2O - 1) * 1000

cPRECT_H218O = cnc.variable("PRECT_H218O", box)
cPRECT_H2O  = cnc.variable("PRECT_H2O", box)
cRp = (cPRECT_H218O / cPRECT_H2O - 1) * 1000

cQFLX_H218O = cnc.variable("QFLX_H218O", box)
cQFLX_H2O = cnc.variable("QFLX_H2O", box)
cRe = (cQFLX_H218O / cQFLX_H2O - 1) * 1000


# 1) Re_bar
# Re_bar = (tRe + cRe) / 2
Re_bar = (tRe + cRe) / 2

# 2) del_E
del_E = tQFLX_H2O - cQFLX_H2O

# 3) E_bar
E_bar = (tQFLX_H2O + cQFLX_H2O) / 2

# 4) del_Re
del_Re = tRe - cRe

# 8) del_Rc
tRc = (tPRECT_H218O - tQFLX_H218O) / (tPRECT_H2O - tQFLX_H2O)
cRc = (cPRECT_H218O - cQFLX_H218O) / (cPRECT_H2O - cQFLX_H2O)
del_Rc = (tRc - cRc) * 1000

# 5) Rc_bar
Rc_bar = (((tRc - 1) * 1000) + ((cRc - 1) * 1000)) / 2 

# 6) del_C
tC = tPRECT_H2O - tQFLX_H2O
cC = cPRECT_H2O - cQFLX_H2O
del_C = tC - cC

# 7) C_bar
C_bar = (tC + cC) / 2

# 10) del_P
del_P = tPRECT_H2O - cPRECT_H2O

# 11) P_bar
P_bar = (tPRECT_H2O + cPRECT_H2O) / 2

# 9) Rp_bar
Rp_bar = (tRp + cRp) / 2 


##################
# Combine it all #
##################

one = Re_bar*del_E
two = E_bar*del_Re 
three = Rc_bar*del_C
four = C_bar*del_Rc 
five = -Rp_bar*del_P
SUM = one + two + three + four + five

print([x/SUM for x in [one, two, three, four, five]])

del_Rp = SUM / P_bar




#-------#



# Done with ratios 
tPRECT_H218O = tnc.variable("PRECT_H218O", box)
tPRECT_H2O  = tnc.variable("PRECT_H2O", box)
tRp = tPRECT_H218O / tPRECT_H2O

tQFLX_H218O = tnc.variable("QFLX_H218O", box)
tQFLX_H2O = tnc.variable("QFLX_H2O", box)
tRe = tQFLX_H218O / tQFLX_H2O

cPRECT_H218O = cnc.variable("PRECT_H218O", box)
cPRECT_H2O  = cnc.variable("PRECT_H2O", box)
cRp = cPRECT_H218O / cPRECT_H2O

cQFLX_H218O = cnc.variable("QFLX_H218O", box)
cQFLX_H2O = cnc.variable("QFLX_H2O", box)
cRe = cQFLX_H218O / cQFLX_H2O


# 1) Re_bar
# Re_bar = (tRe + cRe) / 2
Re_bar = (tRe + cRe) / 2

# 2) del_E
del_E = tQFLX_H2O - cQFLX_H2O

# 3) E_bar
E_bar = (tQFLX_H2O + cQFLX_H2O) / 2

# 4) del_Re
del_Re = tRe - cRe

# 8) del_Rc
tRc = (tPRECT_H218O - tQFLX_H218O) / (tPRECT_H2O - tQFLX_H2O)
cRc = (cPRECT_H218O - cQFLX_H218O) / (cPRECT_H2O - cQFLX_H2O)
del_Rc = tRc - cRc

# 5) Rc_bar
Rc_bar = (tRc + cRc) / 2 

# 6) del_C
tC = tPRECT_H2O - tQFLX_H2O
cC = cPRECT_H2O - cQFLX_H2O
del_C = tC - cC

# 7) C_bar
C_bar = (tC + cC) / 2

# 10) del_P
del_P = tPRECT_H2O - cPRECT_H2O

# 11) P_bar
P_bar = (tPRECT_H2O + cPRECT_H2O) / 2

# 9) Rp_bar
Rp_bar = (tRp + cRp) / 2 


##################
# Combine it all #
##################

one = Re_bar*del_E
two = E_bar*del_Re 
three = Rc_bar*del_C
four = C_bar*del_Rc 
five = -Rp_bar*del_P
SUM = one + two + three + four + five

print([x/SUM for x in [one, two, three, four, five]])

del_Rp = SUM / P_bar