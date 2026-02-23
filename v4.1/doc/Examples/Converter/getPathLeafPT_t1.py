# - getPathLeaf (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

p = Internal.getPathLeaf('CGNSTree/Base/Zone')
test.testO(p, 1)
