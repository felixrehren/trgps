
#
#  transposition groups implementation
#

InstallMethod( DuplicateGroup,
	[IsGroup],
	function( G )
	local H, ktP, x;
	H := GroupByGenerators(GeneratorsOfGroup(G));
	for x in List(KnownAttributesOfObject(G),EvalString)
	do Setter(x)(H,x(G)); od;
	ktP := KnownTruePropertiesOfObject(G);
	for x in List(ktP,EvalString)
	do Setter(x)(H,true); od;
	for x in List(Difference(KnownPropertiesOfObject(G),ktP),EvalString)
	do Setter(x)(H,false); od;
	return H;
	end
);

InstallMethod( TranspositionGroup,
	[IsGroup,IsCollection],
	function( G, D )
	return TrgpNC(G,List(Orbits(G,D),Representative));
	end
	);
	InstallMethod( TrgpNC,
	[IsGroup,IsCollection],
	function( G, D )
	local H;
	H := DuplicateGroup(G);
	ResetFilterObj( H, HasTranspositions );
	ResetFilterObj( H, HasPairs );
	ResetFilterObj( H, HasIncidencePairs );
	ResetFilterObj( H, HasDihedrals );
	ResetFilterObj( H, HasMinimalDihfusController );
	SetTranspositions(H,D);
	ResetFilterObj( H, HasAutomorphismGroup );
	Unbind( H!.AutomorphismGroup ); # forcibly remove
	return H;
	end
	);
InstallMethod( StringTrpoClasses@,
	[IsTrgp],
	function(T)
	return JoinStringsWithSeparator(
		List(Transpositions(T),t->OrbitLength(T,t)),"+");
	end
	);
InstallMethod( PrintString, "for a trgp",
	[IsTrgp],
	function( T )
	local gptxt;
	ResetFilterObj(T,IsTrgp);
	gptxt := PrintString(T);
	SetFilterObj(T,IsTrgp);
	return Concatenation(
		"TrgpNC(\n",
		"\t",gptxt,",\n",
		"\t",PrintString(Transpositions(T)),"\n",
		")" );
	end
	);
InstallMethod( Description, "for a trgp",
	[IsTrgp],
	function( T )
	return Concatenation(
		StructureDescription(T:short),
		" with ",
		StringTrpoClasses@(T),
		" transpos"
	);
	end
	);
	InstallMethod( ViewString, "for a trgp",
	[IsTrgp],
	Description
	);
InstallMethod( \=,
	[IsTrgp,IsTrgp],
	function( A, B )
	ResetFilterObj(A,IsTrgp);
	ResetFilterObj(B,IsTrgp);
	if A <> B then return false;
	else
		SetFilterObj(A,IsTrgp);
		SetFilterObj(B,IsTrgp);
		return ForAll(Transpositions(A),t->
			ForAny(Transpositions(B),s->IsConjugate(A,t,s)) );
	fi;
	end
);

InstallMethod( CartWoDiag@,
	[IsCollection,IsCollection],
	function( A, B )
	local r, i, j;
	r := [];
	for i in [1..Length(A)] do for j in [1..Length(B)] do
		if A[i] <> B[j] then
			Add(r,Set([A[i],B[j]])); fi; od; od;
	return r;
	end
	);
InstallMethod( Pairs,
	[IsTrgp],
	function( T )
  local tt;
	tt := Union(List(Transpositions(T),t->AsList(t^T)));
	return List(
		OrbitsDomain(
			T,
			CartWoDiag@( tt,tt ),
			OnSets),
		Representative);
	return Pairs(T);
	end
	);
InstallMethod( TranspositionDegree,
	[IsTrgp],
	T -> Maximum(List(Pairs(T),p->Order(Product(p))))
);

InstallMethod( IsInvolution,
	[IsMultiplicativeElementWithOne],
	g -> IsOne(g^2)
	and not IsOne(g)
	);
InstallMethod( IsOrderIn@,
	[IsMultiplicativeElementWithOne,IsList],
	function( g, N )
	return ForAny( N, n-> IsOne(g^n) ); end
	);
InstallMethod( CanBeTrgp,
	[IsGroup,IsList],
	function( G, N )
  local C, cc, mat, PairsOrderIn, plant, seeds, L, i, j, l;
	C := CommutatorFactorGroup(G);
	if not ((Size(C) mod 2 = 0) and IsElementaryAbelian(C))
		then return false; fi;
	cc := Filtered( ConjugacyClasses( G ),
		c -> IsInvolution(Representative(c)) );
	if G <> Subgroup(G,Concatenation(List(cc,AsList))) then return false; fi;
	mat := List([1..Length(cc)],i->[]);
	PairsOrderIn := function( S,T, N )
		return ForAll( Cartesian(S,T),st->IsOrderIn@(Product(st),N) ); end;
	for i in [1..Length(cc)] do for j in [i..Length(cc)] do
		mat[i][j] := PairsOrderIn(AsList(cc[i]),AsList(cc[j]),N);
		mat[j][i] := mat[i][j]; od; od;
	plant := Filtered([1..Length(cc)],i->ForAll(mat[i],IdFunc));
	if G = Subgroup(G,Concatenation(List(cc{plant},AsList)))
		then return true; fi; # shortcut?
	seeds := Filtered([1..Length(cc)],i->mat[i][i] and not (i in plant));
	L := List(seeds,i->Union(plant,[i]));
	if ForAny(L,l->G = Subgroup(G,Concatenation(List(cc{l},AsList))))
		then return true; fi;
	for i in seeds do
		for l in L do
			if i > Maximum(l)
			and ForAll(l,j->mat[j][i])
			then Add(L,Concatenation(l,[i]));
				if G = Subgroup(G,Concatenation(List(cc{Concatenation(l,[i])},AsList)))
				then return true; fi;
			fi;
		od;
	od;
	return false;
	end
	);
	InstallMethod( CanBeTrgp,
	[IsGroup,IsPosInt],
	function( G, N )
	return CanBeTrgp( G, [2..N] ); end
	);
InstallMethod( GroupToTrgps,
	[IsGroup,IsList],
	function( G, N )
  local cc, pairs, gph, okCC, cliques, newc, S, OnCC, OnCliq, Transpositions, i, j, minlOnly, k, c;
	cc := Filtered(
		ConjugacyClasses(G),
		x -> IsInvolution(Representative(x)) );
	pairs := List([1..Length(cc)],i->[]);
	gph := List([1..Length(cc)],i->[]);
	for i in [1..Length(cc)] do for j in [1..i] do
			pairs[i][j] := List(
				OrbitsDomain(
					G,
					CartWoDiag@(AsList(cc[i]),AsList(cc[j])),
					OnSets ),
				o -> o[1] );
			pairs[j][i] := pairs[i][j];
			gph[i][j] := ForAll( pairs[i][j],
				p -> IsOrderIn@(Product(p),N) );
			gph[j][i] := gph[i][j];
	od; od;
	okCC:= List(
		Filtered(
			[1..Length(cc)],
			i -> gph[i][i] ),
		i -> rec(
			class:=i,
			size:=Size(Subgroup(G,AsList(cc[i]))) )
		);
	cliques := List(
		okCC,
		k -> rec(
			clique:=[k.class],
			size:=k.size )
		);
	minlOnly := (ValueOption("minimal") = true);
	for k in okCC do
		for c in cliques do
			if not (minlOnly and (c.size=Size(G)))
				and k.class > Maximum(c.clique)
				and ForAll(c.clique,i->gph[k.class][i])
			then
				newc := Concatenation(c.clique,[k.class]);
				if c.size = Size(G) or k.size = Size(G) then S := Size(G);
				else S := Size(Subgroup(G,Concatenation(List(cc{newc},AsList)))); fi;
				if (not minlOnly or S > c.size)
				then Add(cliques, rec(clique:=newc,size:=S)); fi;
			fi;
		od;
	od;
	cliques := List(Filtered(cliques,c->c.size=Size(G)),c->c.clique);
	if minlOnly then for c in cliques
	do cliques := cliques{Filtered([1..Length(cliques)],
			i->c=cliques[i] or not IsSubset(cliques[i],c))};
	od; fi;
	OnCC := function( i, a )
		return Position(cc,(Representative(cc[i])^a)^G); end;
	OnCliq := function( c, a )
		local d;
		d := List(c,i -> OnCC(i,a));
		Sort(d);
		return d; end;
	cliques := List(
		OrbitsDomain(
			AutomorphismGroup(G),
			cliques,
			OnCliq),
		x -> x[1] );
	return List(
		cliques,
		function(c) local T;
		T := TranspositionGroup( G, List(cc{c},Representative) );
		SetPairs(T,Concatenation(Concatenation(
			List([1..Length(c)],i->List([1..i],
			j->pairs[c[i]][c[j]] )) )));
		return T;
		end );
	end
	);
	InstallMethod( GroupToTrgps,
	[IsGroup,IsPosInt],
	function( G, N )
	return GroupToTrgps( G, [2..N] );
	end
	);
	InstallMethod( GroupToTrgps,
	[IsGroup],
	function( G )
	return GroupToTrgps( G, [2..Size(G)] );
	end
	);
InstallMethod( IsMinimalTrgp,
	[IsTrgp],
	function(T)
	local cc;
	if Length(Transpositions(T)) = 1 then return true;
	cc := List(Transpositions(T),t->t^T);
	else return ForAll(
		Combinations(cc,Length(cc)-1),
		c->Group(Union(c)) <> Group(Union(cc))
	); fi; end
	);
InstallMethod( OnTrgps,
	[IsTrgp,IsMultiplicativeElement],
	function( om, g )
	return TranspositionGroup( om^g, List(Transpositions(om),x->x^g) );
	end
	);
InstallMethod( AutomorphismGroup,
	"for transposition groups",
	[IsTrgp],
	SUM_FLAGS,
	function(T)
	local tt, A;
	tt := T!.Transpositions;
	Unbind(T!.Transpositions);
	ResetFilterObj(T,IsTrgp);
	A := AutomorphismGroup(T);
	ResetFilterObj(T,HasAutomorphismGroup);
	SetTranspositions(T,tt);
	return Subgroup(A, Filtered(A, a -> ForAll(Transpositions(T),
		t->ForAny(Transpositions(T),s->IsConjugate(T,t^a,s)) )));
	end
);

InstallMethod( TrgpSearch,
	[IsPosInt,IsPosInt],
	function( s, N )
  local thresh, gg, nn, time, G, n;
	thresh := 20;
	gg := [];
	nn := NumberSmallGroups( s );
	time := Runtime();
	if nn > thresh then
		Print(nn, " groups total\n");
		for n in [1..nn] do
			if n mod thresh = 0 then
				Print("\tchecked ",n," groups (in ",time," ms)\n");
				time := Runtime(); fi;
			if n mod 1000 = 0 then
				Print("#searching through the ",nn," gps",
				" of size ",s," for ",N,"-trgps\n"); fi;
			G := SmallGroup( s,n );
			if CanBeTrgp( G, N ) then
				Add(gg,[s,n]); fi;
		od;
	else
		gg := IdsOfAllSmallGroups(Size,s,G -> CanBeTrgp(G,N),true); fi;
	return gg;
	end
	);
InstallMethod( TrgpQuest,
	[IsList,IsPosInt,IsString],
	function( R, N, file )
  local time, gg, r;
	for r in R do
		if r mod 2 <> 0 then continue; fi;
		time := Runtime();
		Print(r,": ");
		gg := TrgpSearch( r, N );
		Print(Size(gg)," trgps / ",NumberSmallGroups(r)," gps ",
			"(in ",time," ms)\n");
		if not IsEmpty(gg) then
			AppendTo(file,JoinStringsWithSeparator(List(gg,PrintString),", "),",\n"); fi;
		UnloadSmallGroupsData();
	od;
	end
	);
InstallMethod( ViaAtlas@,
	[IsString,IsList],
	function( GName, ccNames )
  local G, info, gens, cclspr, rr, filter;
	LoadPackage("atlasrep", false);
	G := AtlasGroup( GName );
	info := OneAtlasGeneratingSetInfo(GName);
	if info <> fail then
		gens := AtlasGenerators( info.identifier );
		cclspr := AtlasProgram( GName, gens.standardization, "classes" );
		if cclspr = fail then
			if IsBound( ccNames )
			then return fail;
			else
				rr := List(
					Filtered( ConjugacyClasses( G ),
						c->IsInvolution(Representative(c)) ),
					Representative);
			fi;
		else
			if not IsEmpty( ccNames ) then
				filter := i -> cclspr.outputs[i] in ccNames;
			else
				filter := i -> cclspr.outputs[i][1] = '2';
			fi;
			cclspr := RestrictOutputsOfSLP( cclspr.program,
								Filtered( [1..Length( cclspr.outputs )], filter ) );
			rr := ResultOfStraightLineProgram( cclspr, gens.generators );
		fi;
	fi;
	return TranspositionGroup( G, rr );
	end
);

InstallMethod( IsIsomOfTrgp,
	[IsTrgp,IsTrgp,IsMapping],
	function(R,S,h)
	if R <> Source(h) or S <> Image(h)
	or Size(R) <> Size(S)
	or Length(Transpositions(R)) <> Length(Transpositions(S))
	then return false;
	else return 
		ForAll( Transpositions(R), t -> ForAny( Transpositions(S), u ->
			IsConjugate(S,Image(h,t),u)) );
	fi;
	end
	);
InstallMethod( AreIsomorphicTrgp,
	[IsTrgp,IsTrgp],
	function( R, S )
	if Size(R) <> Size(S)
	or Length(Transpositions(R)) <> Length(Transpositions(S))
	then return false; fi;
	return ForAny( AllHomomorphismClasses( R, S ),
		h -> IsIsomOfTrgp(R,S,h) );
	end
	);
InstallMethod( IsomorphismClassesTrgps,
	[IsTrgp,IsTrgp],
	function( R, S )
	if Size(R) <> Size(S)
	or Length(Transpositions(R)) <> Length(Transpositions(S))
	then return []; fi;
	return Filtered( AllHomomorphismClasses( R, S ),
		h -> IsIsomOfTrgp(R,S,h) );
	end
	);
InstallMethod( EmbeddingsClassesTrgps,
	[IsTrgp,IsTrgp],
	function( R, S )
	return Filtered(
		AllHomomorphismClasses( R, S ),
		h -> IsInjective(h)
			and IsIsomOfTrgp(R,TrgpNC(Image(h),Image(h,Transpositions(R))),h)
		);
	end
	);
InstallMethod( MaximalSubgroupsTrgp,
	[IsTrgp],
	function( R )
  local cc, ms, i, j;
	cc := Concatenation(List(Transpositions(R),t->AsList(t^R)));
	ms := List(
		MaximalSubgroupClassReps(R),
		m -> Subgroup(R,Intersection(m,cc)) );
	for i in [Length(ms),Length(ms)-1..1] do
		if Size(ms[i]) < 3 then Remove(ms,i); continue; fi;
		for j in [1..i-1] do
			if IsConjugate(R,ms[i],ms[j]) then
				Remove(ms,i); break; fi;
	od; od;
	return List(
		ms,
		sgp -> TranspositionGroup(
			sgp,
			Intersection(sgp,cc) )
		);
	end
	);
InstallMethod( MaximalSubtrposTrgp,
	[IsTrgp],
	R ->
	List(
		Filtered(
			Combinations(
				Transpositions(R),
				Length(Transpositions(R))-1
			),
			c -> R = Subgroup(R,Concatenation(List(c,k->AsList(k^R))))
		),
		c -> TranspositionGroup(R,c)
	)
	);
InstallMethod( Subtrgp,
	[IsTrgp, IsGroup],
	function( T, H )
  local D;
	D := Intersection(Concatenation(List(Transpositions(T),t->t^T)),H);
	if H <> Subgroup(T,D) then return fail; fi;
	return TrgpNC(H,List(OrbitsDomain(H,D),Representative));
	end
	);
InstallMethod( ImagesSet,"for a transposition group",
	[IsMapping,IsTrgp],SUM_FLAGS,
	function( f, S )
  return TranspositionGroup(
		Image(f),
		Image(f,Transpositions(S)) );
	end
	);
InstallMethod( AsSmallPermTrgp,
	[IsTrgp],
	function(T)
	local S;
	if IsPermGroup(T)
	then 
		S := Images(SmallerDegreePermutationRepresentation(T),T);
		S!.Pairs := List(S!.Pairs,Set);
		return S;
	else return AsSmallPermTrgp(Images(IsomorphismPermGroup(T),T)); fi;
	end
);

InstallMethod( IncidencePairs,
	[IsTrgp],
	function( T )
	local L, l, silhm, GivesSubAlg, i, j, ArrangePartsPerm, phi;
	L := List(Pairs(T),x->Order(x[1]*x[2]));
	l := Length(L);
	silhm := IdentityMat(l);
	GivesSubAlg := function(G,p1,p2)
		local g1,g2,s1,s2;
		g1 := Group(p1);
		s1 := Union(List(p1,p->Orbit(g1,p)));
		g2 := Group(p2);
		s2 := Union(List(p2,p->Orbit(g2,p)));
		return ForAny(G,g->IsSubset(s2,OnSets(s1,g)));
	end;
	for i in [1..l] do for j in [1..l] do
		if L[i] <= L[j] and IsInt(L[j]/L[i])
		and GivesSubAlg(T,Pairs(T)[i],Pairs(T)[j])
		then silhm[i][j] := 1; fi;
	od; od;

	ArrangePartsPerm := function(M,orders)
		local IncPart, part, U, sg, phi, PF, PL, V, P;
		IncPart := function(M)
			local part, L, Comp, C, i, c;
			part := [];
			L := List([1..Length(M)],i->true);
			Comp := function(i)
				local C, j,x, X, Y, Z;
				C := [i];
				for j in C do
					X := Concatenation([1..j-1],[j+1..Length(M)]);
					Y := Filtered(X, x-> IsOne(M[j][x]));
					Z := Filtered(X, x-> IsOne(M[x][j]));
					for x in Union(Y,Z) do
						if not x in C then
							Add(C,x);
						fi;
					od;
				od;
				return C;
			end;
			for i in [1..Length(L)] do
				if L[i] then
					C := Comp(i);
					for c in C do
						L[c] := false; od;
					Add(part, C);
				fi;
			od;
			return part;
		end;
		part := IncPart(M);
		U := Concatenation(part);
		sg := PermListList(U,[1..Length(M)]);
		M := Permuted(M,sg);
		Apply(M, m->Permuted(m,sg));
		phi := ();
		for P in part do
			PF := Position(U, P[1]);
			PL := Position(U, P[Length(P)]);
			V := List(U, function(u)
						if   Position(U,u) < PF then return Position(U,u)-100;
						elif Position(U,u) > PL then return Position(U,u)+100;
						else return orders[u]; fi; end);
			phi := phi*SortingPerm(V);
		od;
		return sg*phi;
	end;
	phi := ArrangePartsPerm(silhm,L);
	T!.Pairs := Permuted(Pairs(T),phi);
	silhm := Permuted(silhm,phi);
	Apply(silhm, s->Permuted(s,phi));
	return silhm;
	end
);

InstallMethod( Dihedrals,
	[IsTrgp],
	T -> 
	Set(List( Pairs(T),p -> Group(p)^T ))
	);
InstallMethod( ControlsFusion,
	[IsCollection,IsGroup],
	function( Dl, H )
	return ForAll( Dl, dl -> Length(OrbitsDomain(H,dl))=1 );
	end
	);
InstallMethod( ControlsDihedralFusion,
	[IsTrgp,IsGroup],
	function( T, H )
	return ControlsFusion( Dihedrals(T), H );
	end
	);
InstallMethod( MinimalFusionController,
	[IsCollection,IsGroup],
	function( Dl, H )
	local res, M;
	res := [];
	for M in MaximalSubgroupClassReps(H) do
		if ControlsFusion( Dl, M )
		then Append(res,MinimalFusionController( Dl, M )); fi;
	od;
	if IsEmpty(res) then return [H];
	else return res; fi;
	end
	);
	InstallMethod( MinimalDihfusController,
	[IsTrgp],
	T -> MinimalFusionController(Dihedrals(T),T)
);
