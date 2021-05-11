#ifndef _MYMAP_H
#define _MYMAP_H
#include <map>
#include <set>
#include <string>
/// \cond
namespace ZMesh {
/// \endcond
class myComp
{
	public:
		bool operator()( std::string s1, std::string s2)
		{
			return s1<s2;
		}
};

class patchType;
class ZPatch;
class ZPatch2;
class ZElem;
class ZNode;
class ZEntityLessThen;
class nodeField;
class Eqn;

typedef std::map<std::string, patchType, myComp> str2PatchMap; 
typedef std::map<std::string, ZPatch, myComp> str2ZPatchMap; 
typedef std::map<std::string, ZPatch2, myComp> str2ZPatch2Map; 
typedef std::map<std::string, int, myComp> str2IntMap;
typedef std::map<int, std::string> int2StrMap;
typedef std::set<ZElem*, ZEntityLessThen> elemSet;
typedef std::set<ZNode*, ZEntityLessThen> nodeSet;
typedef std::map<std::string, nodeField*, myComp> str2FieldPointer;
typedef std::map<std::string, Eqn*, myComp> str2EqnPointer;

} // end namespace ZMesh
#endif
