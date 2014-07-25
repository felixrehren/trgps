
#
#  transposition groups declarations
#

DeclareOperation("DuplicateGroup",[IsGroup]);

DeclareAttribute("Transpositions",IsGroup);
DeclareSynonym("IsTrgp",HasTranspositions);
DeclareOperation("Description",[IsTrgp]);

DeclareOperation("TranspositionGroup",[IsGroup,IsCollection]);
DeclareSynonym("Trgp", TranspositionGroup);
DeclareOperation("TrgpNC",[IsGroup,IsCollection]);
DeclareOperation("StringTrpoClasses@",[IsTrgp]);

DeclareOperation("CartWoDiag@",[IsCollection,IsCollection]);
DeclareAttribute("Pairs",IsTrgp);
DeclareAttribute("TranspositionDegree",IsTrgp);
DeclareAttribute("Regularity",IsTrgp);

DeclareOperation("IsInvolution",[IsMultiplicativeElementWithOne]);
DeclareOperation("IsOrderIn@",[IsMultiplicativeElementWithOne,IsList]);
DeclareOperation("CanBeTrgp",[IsGroup]);
DeclareOperation("CanBeTrgp",[IsGroup,IsList]);
DeclareOperation("CanBeTrgp",[IsGroup,IsPosInt]);
DeclareOperation("GroupToTrgps",[IsGroup]);
DeclareOperation("GroupToTrgps",[IsGroup,IsList]);
DeclareOperation("GroupToTrgps",[IsGroup,IsPosInt]);
DeclareOperation("IsMinimalTrgp",[IsTrgp]);
DeclareOperation("OnTrgps",[IsTrgp,IsMultiplicativeElement]);

DeclareOperation("TrgpSearch",[IsPosInt,IsPosInt]);
DeclareOperation("TrgpQuest",[IsList,IsPosInt,IsString]);
DeclareOperation("ViaAtlas@",[IsString,IsList]);

DeclareOperation("IsIsomOfTrgp",[IsTrgp,IsTrgp,IsMapping]);
DeclareOperation("AreIsomorphicTrgp",[IsTrgp,IsTrgp]);
DeclareOperation("IsomorphismClassesTrgps",[IsTrgp,IsTrgp]);
DeclareOperation("EmbeddingsClassesTrgps",[IsTrgp,IsTrgp]);
DeclareOperation("MaximalSubgroupsTrgp",[IsTrgp]);
DeclareOperation("MaximalSubtrposTrgp",[IsTrgp]);
DeclareOperation("Subtrgp",[IsTrgp,IsGroup]);
DeclareOperation("AsSmallPermTrgp",[IsTrgp]);

DeclareAttribute("IncidencePairs",IsTrgp);
DeclareOperation("PairsGraph",[IsTrgp]);
DeclareOperation("PairsGraph",[IsMatrix,IsList]);

DeclareOperation("NoncommGraph",[IsTrgp]);
DeclareOperation("NoncommGraph",[IsList]);

DeclareOperation("CoxeterGroup",[IsMatrix]);
DeclareOperation("CoxeterGroup",[IsRootSystem]);
DeclareOperation("CoxeterGroup",[IsString,IsPosInt]);

DeclareAttribute("Dihedrals",IsGroup);
DeclareOperation("ControlsFusion",[IsCollection,IsGroup]);
DeclareOperation("ControlsDihedralFusion",[IsTrgp,IsGroup]);
DeclareOperation("MinimalFusionController",[IsCollection,IsGroup]);
DeclareAttribute("MinimalDihfusController",IsTrgp);
