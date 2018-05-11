PRECT=0
PRECC=0
PRECL=0
PRECSC=0
PRECSL=0
PRECST=0
QFLX=0
CLDHGH=0
PS=0
PSL=0
PRECT_d18O=0
PRECT_dD=0
PRECT_dxs=0
Q850=0
Q500=0
Q200=0
V850=0
V500=0
V200=0
VT850=0
VT500=0
VT200=0
VQ850=0
VQ500=0
VQ200=0
U850=0
U500=0
U200=0
UT850=0
UT500=0
UT200=0
UQ850=0
UQ500=0
UQ200=0
T850=0
T500=0
T200=0
OMEGA850=0
OMEGA500=0
OMEGA200=0
Z3850=0
Z3500=0
Z3200=0
dDV850=0
dDV500=0
dDV200=0
d18OV850=0
d18OV500=1
d18OV200=0
dxsV850=0
dxsV500=0
dxsV200=0












#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create variable string list to loop through
var2dlist = []
if PRECT:
	var2dlist.append("PRECT")
if PRECC:
	var2dlist.append("PRECC")
if PRECL:
	var2dlist.append("PRECL")
if PRECSC:
	var2dlist.append("PRECSC")
if PRECSL:
	var2dlist.append("PRECSL")
if PRECST:
	var2dlist.append("PRECST")
if QFLX:
	var2dlist.append("QFLX")
if CLDHGH:
	var2dlist.append("CLDHGH")
if PS:
	var2dlist.append("PS")
if PSL:
	var2dlist.append("PSL")
if PRECT_d18O:
	var2dlist.append("PRECT_d18O")
if PRECT_dD:
	var2dlist.append("PRECT_dD")
if PRECT_dxs:
	var2dlist.append("PRECT_dxs")

var3dlist = []
if Q850:
	var3dlist.append("Q850")
if Q500:
	var3dlist.append("Q500")
if Q200:
	var3dlist.append("Q200")
if V850:
	var3dlist.append("V850")
if V500:
	var3dlist.append("V500")
if V200:
	var3dlist.append("V200")
if VT850:
	var3dlist.append("VT850")
if VT500:
	var3dlist.append("VT500")
if VT200:
	var3dlist.append("VT200")
if VQ850:
	var3dlist.append("VQ850")
if VQ500:
	var3dlist.append("VQ500")
if VQ200:
	var3dlist.append("VQ200")
if U850:
	var3dlist.append("U850")
if U500:
	var3dlist.append("U500")
if U200:
	var3dlist.append("U200")
if UT850:
	var3dlist.append("UT850")
if UT500:
	var3dlist.append("UT500")
if UT200:
	var3dlist.append("UT200")
if UQ850:
	var3dlist.append("UQ850")
if UQ500:
	var3dlist.append("UQ500")
if UQ200:
	var3dlist.append("UQ200")
if T850:
	var3dlist.append("T850")
if T500:
	var3dlist.append("T500")
if T200:
	var3dlist.append("T200")
if OMEGA850:
	var3dlist.append("OMEGA850")
if OMEGA500:
	var3dlist.append("OMEGA500")
if OMEGA200:
	var3dlist.append("OMEGA200")
if Z3850:
	var3dlist.append("Z3850")
if Z3500:
	var3dlist.append("Z3500")
if Z3200:
	var3dlist.append("Z3200")
if dDV850:
	var3dlist.append("dDV850")
if dDV500:
	var3dlist.append("dDV500")
if dDV200:
	var3dlist.append("dDV200")
if d18OV850:
	var3dlist.append("d18OV850")
if d18OV500:
	var3dlist.append("d18OV500")
if d18OV200:
	var3dlist.append("d18OV200")
if dxsV850:
	var3dlist.append("dxsV850")
if dxsV500:
	var3dlist.append("dxsV500")
if dxsV200:
	var3dlist.append("dxsV200")