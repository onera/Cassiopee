import os
paths1 = '../ThirdParty/paradigma/paradigm/src'
paths2 = 'XCore/paradigma' # local link
if not os.path.exists(paths2) and os.path.exists(paths1):
    os.symlink('../'+paths1, paths2)  
