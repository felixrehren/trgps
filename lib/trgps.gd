
#
#  transposition groups declarations
#

FamilyTrgp@ := NewFamily("TranspositionGroupsFamily");
DeclareCategory("IsTrgp",IsObject);
DeclareRepresentation("IsTrgpStdRep",
	IsComponentObjectRep and IsAttributeStoringRep,[]);
TypeTrgp@ := NewType(FamilyTrgp@,IsTrgp and IsTrgpStdRep);

DeclareOperation("TranspositionGroup",[IsGroup,IsCollection]);
DeclareOperation("Trgp",[IsGroup,IsCollection]);
DeclareOperation("Transpositions",[IsTrgp]);
DeclareOperation("GroupX",[IsTrgp]);
DeclareOperation("StringTrpoClasses@",[IsTrgp]);

DeclareOperation("CartWoDiag@",[IsCollection,IsCollection]);
DeclareAttribute("Pairs",IsTrgp);
DeclareOperation("TranspositionDegree",[IsTrgp]);

DeclareOperation("IsInvolution",[IsMultiplicativeElementWithOne]);
DeclareOperation("IsOrderIn@",[IsMultiplicativeElementWithOne,IsList]);
DeclareOperation("CanBeTrgp",[IsGroup,IsList]);
DeclareOperation("CanBeTrgp",[IsGroup,IsPosInt]);
DeclareOperation("GroupToTrgps",[IsGroup,IsList,IsBool]);
DeclareOperation("GroupToTrgps",[IsGroup,IsList]);
DeclareOperation("GroupToTrgps",[IsGroup,IsPosInt,IsBool]);
DeclareOperation("GroupToTrgps",[IsGroup,IsPosInt]);
DeclareOperation("IsMinimalTrgp",[IsTrgp]);
DeclareOperation("OnTrgps",[IsTrgp,IsMultiplicativeElement]);
DeclareAttribute("AutomorphismGroup",IsTrgp);

DeclareOperation("TrgpSearch",[IsPosInt,IsPosInt]);
DeclareOperation("TrgpQuest",[IsList,IsPosInt,IsString]);
DeclareOperation("ViaAtlas@",[IsString,IsList]);

DeclareOperation("IsomorphismClassesTrgps",[IsTrgp,IsTrgp]);
DeclareOperation("EmbeddingsClassesTrgps",[IsTrgp,IsTrgp]);
DeclareOperation("MaximalSubgroupsTrgp",[IsTrgp]);
DeclareOperation("MaximalSubtrposTrgp",[IsTrgp]);
DeclareOperation("Subtrgp",[IsTrgp,IsGroup]);
DeclareOperation("ImageTrgp",[IsTrgp,IsMapping]);
DeclareOperation("TrgpSmallerDegreeRep",[IsTrgp]);

DeclareAttribute("IncidencePairs",IsTrgp);

DeclareAttribute("Dihedrals",IsTrgp);
DeclareOperation("ControlsDihedralFusion",[IsCollection,IsGroup]);
DeclareOperation("ControlsDihedralFusion",[IsTrgp,IsGroup]);
DeclareOperation("MinimalDihfusController",[IsCollection,IsGroup]);
DeclareOperation("MinimalDihfusController",[IsTrgp]);
