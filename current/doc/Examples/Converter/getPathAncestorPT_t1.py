# - getPathAncestor (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

p = Internal.getPathAncestor('CGNSTree/Base/Zone', level=1)
test.testO(p, 1)
