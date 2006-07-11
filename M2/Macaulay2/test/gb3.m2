R = QQ[x,y,z,MonomialOrder=>Lex]
gb ideal {x^5+y^4+z^3-1, x^3+y^2+z^2-1}


S = QQ[a..d]
I = ideal(a^2-1, a^3+2)
T = S/I
assert( T == 0 )

-- May 23, 2006:
R = ZZ[s,t,x,y,z, MonomialOrder=>{2,3}];
I = ideal(x-s^3-s*t-1, y-t^3-3*t^2-t, z-s*t^3)
leadTerm gens gb I  --crashes on my macintel.

gens(g = gb(f = random (ZZ^6, ZZ^6), ChangeMatrix => true))
assert( gens g == f * getChangeMatrix g )
assert( leadTerm gens g == transpose leadTerm gens g)

-- uniqueness of gb's:
assert( gens gb matrix {{-1}} == gens gb matrix {{1}} )

-- gb over ZZ
f = matrix {
     {13650502, 198662, -1514226}, {-528389638951, -7688266050, 58613349522}, {1819050, 26473, -201784},
     {-34721130542, -505205335, 3851555009}, {13943863165, 202888407, -1546768644}, {112371429966, 1635046125, -12465168534}}
h = gens gb f
assert(gens gb h == h)
assert( h == matrix {
	  {-304371114, 129997874, -86994465284}, {-224654892, 95950799, 337132988376}, {-38621018, 16495160, -10628621651},
	  {51393451, -21950302, -74040381498}, {0, 1, -18062908713}, {0, 0, 3}})
assert isIsomorphism map(coker f, coker h, id_(target f))
G = gb(f,ChangeMatrix => true)
c = getChangeMatrix G
assert( h == gens G )
assert isIsomorphism map(image f, image h, c)

end
-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/test gb3.out"
-- End:
