# - newFamilyName (pyTree) -
import Converter.Internal as Internal

# Used to tag a zone for ex. The value refers to the name of a Family_t node
n = Internal.newFamilyName(name='FamilyName', value='Wing'); Internal.printTree(n)
#>> ['FamilyName',array('Wing',dtype='|S1'),[0 son],'FamilyName_t']
