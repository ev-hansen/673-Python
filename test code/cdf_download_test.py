from pyspedas.rbsp import hope 

probe = 'a'
yr_str = "2014"  # year
m_str = "06"   # month
d_str = "20"   # day
datatype = "pitchangle"
datatype_abbr = "pa"
level = "l3"
instrument = "ect"
component = "hope"
rel = "rel04"


hope(trange=[f"{yr_str}-{m_str}-{d_str}", f"{yr_str}-{m_str}-{d_str}"], 
     probe=probe, datatype=datatype, downloadonly=True)
