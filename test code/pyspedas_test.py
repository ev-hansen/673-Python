import pprint
pp = pprint.PrettyPrinter(indent=4)

import pyspedas

yr_str = "2015"
m_str = "11"
d_str = "26"
probe = 'a'
pyspedas.rbsp.hope(trange=[f"{yr_str}-{m_str}-{d_str}", 
                           f"{yr_str}-{m_str}-{d_str}"], probe=probe)  
print("\n\n")
pp.pprint(globals())
print("\n\n")
pp.pprint(locals())
