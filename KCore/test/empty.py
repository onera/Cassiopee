# - empty (numpy like) -
import KCore

# Return un numpy aligne (fortran, double)
a = KCore.empty( (1200,3), 64 )
print a.shape
print a.flags
