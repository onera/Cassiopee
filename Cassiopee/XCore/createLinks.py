import os
paths1 = '../ThirdParty/paradigma/paradigm/src'
paths2 = 'XCore/paradigma' # local link
if not os.path.exists(paths2) and os.path.exists(paths1):
    # symlink fails on windows
    try:
        os.symlink('../'+paths1, paths2)
    except:
        import shutil
        shutil.copytree(paths1, paths2)
