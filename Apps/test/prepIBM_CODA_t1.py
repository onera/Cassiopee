# - Fast.IBM -
import Apps.Coda.ToolboxIBM_CODA as AppIBM
import KCore.test as test
LOCAL = test.getLocal()


# Prepare
AppIBM.prepare('naca1DRANS.cgns', t_out=LOCAL+'/t.cgns',check=True)
