from Qpyl.core.qparameter import QPrm, QPrmError


parameters = QPrm('oplsaa')
parameters.read_prm('qoplsaa.prm')
parameters.read_prm('PAF.prm')
parameters.read_prm('IND.prm')
parameters.read_prm('IMI.prm')
parameters.read_prm('WAT.prm')
parameters.read_prm('PRD.prm')

with open('qoplsaa_all.prm', 'w') as f:
    f.write(parameters.get_string())
    
