// 
// Copyright (c) 2006, Benjamin Kaufmann
// 
// This file is part of Clasp. See http://www.cs.uni-potsdam.de/clasp/ 
// 
// Clasp is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Clasp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Clasp; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>
#include <algorithm>
#include <clasp/weight_constraint.h>
#include <clasp/solver.h>
using namespace std;

namespace Clasp { namespace Test {
class CardTest : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(CardTest);
	CPPUNIT_TEST(testAssertTriviallySat);
	CPPUNIT_TEST(testAssertTriviallyUnSat);
	CPPUNIT_TEST(testAssertNotSoTriviallySat);
	CPPUNIT_TEST(testAssertNotSoTriviallyUnSat);

	CPPUNIT_TEST(testTrivialBackpropTrue);
	CPPUNIT_TEST(testTrivialBackpropFalse);
	CPPUNIT_TEST(testTrivialBackpropFalseWeight);
	
	CPPUNIT_TEST(testForwardTrue);
	CPPUNIT_TEST(testForwardFalse);
	CPPUNIT_TEST(testBackwardTrue);
	CPPUNIT_TEST(testBackwardFalse);

	CPPUNIT_TEST(testForwardTrueConflict);
	CPPUNIT_TEST(testForwardFalseConflict);
	CPPUNIT_TEST(testBackwardTrueConflict);
	CPPUNIT_TEST(testBackwardFalseConflict);

	CPPUNIT_TEST(testReasonBug);
	CPPUNIT_TEST(testWeightReasonAfterBackprop);
	CPPUNIT_TEST(testOrderBug);
	CPPUNIT_TEST(testBackwardAfterForward);
	
	CPPUNIT_TEST(testSimplify);
	CPPUNIT_TEST(testSimplifyCardinality);
	CPPUNIT_TEST(testSimplifyWeight);

	CPPUNIT_TEST(testAssertWeightTriviallySat);
	CPPUNIT_TEST(testAssertWeightTriviallyUnSat);
	CPPUNIT_TEST(testAssertWeightNotSoTriviallySat);
	CPPUNIT_TEST(testAssertWeightNotSoTriviallyUnSat);
	CPPUNIT_TEST(testWeightForwardTrue);
	CPPUNIT_TEST(testWeightForwardFalse);
	CPPUNIT_TEST(testWeightBackwardTrue);
	CPPUNIT_TEST(testWeightBackwardFalse);
	CPPUNIT_TEST(testWeightConflict);

	CPPUNIT_TEST(testCloneWeight);
	CPPUNIT_TEST(testCloneWeightShared);
	CPPUNIT_TEST_SUITE_END(); 
public:
	CardTest() {
		body  = posLit(ctx.addVar(Var_t::body_var));
		a     = posLit(ctx.addVar(Var_t::atom_var));
		b     = posLit(ctx.addVar(Var_t::atom_var));
		c     = posLit(ctx.addVar(Var_t::atom_var));
		d     = posLit(ctx.addVar(Var_t::atom_var));
		e     = posLit(ctx.addVar(Var_t::atom_var));
		f     = posLit(ctx.addVar(Var_t::atom_var));
		ctx.startAddConstraints(10);
	}
	Solver& solver() { return *ctx.master(); }
	void testAssertTriviallySat() {
		LitVec lits;
		lits.push_back(body);
		lits.push_back(a);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 0));
		CPPUNIT_ASSERT(solver().isTrue(body));
		CPPUNIT_ASSERT(ctx.numConstraints() == 0);
	}
	void testAssertTriviallyUnSat() {
		LitVec lits;
		lits.push_back(body);
		lits.push_back(a);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 2));
		CPPUNIT_ASSERT(solver().isFalse(body));
		CPPUNIT_ASSERT(ctx.numConstraints() == 0);
	}

	void testAssertNotSoTriviallySat() {
		LitVec lits = makeLits();
		solver().force(lits[1], 0);
		solver().force(lits[2], 0);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 2));
		CPPUNIT_ASSERT(solver().isTrue(body));
		CPPUNIT_ASSERT(ctx.numConstraints() == 0);
	}

	void testAssertNotSoTriviallyUnSat() {
		LitVec lits = makeLits();
		solver().force(~lits[1], 0);
		solver().force(~lits[3], 0);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 3));
		CPPUNIT_ASSERT(solver().isTrue(~body));
		CPPUNIT_ASSERT(ctx.numConstraints() == 0);
	}
	
	void testTrivialBackpropTrue() {
		WeightLitVec lits = makeWeightLits();
		solver().force(body, 0);
		CPPUNIT_ASSERT_EQUAL(true,  WeightConstraint::newWeightConstraint(ctx, body, lits, 7));
		CPPUNIT_ASSERT(ctx.numConstraints() == 1);
		CPPUNIT_ASSERT(solver().isTrue(a));
		CPPUNIT_ASSERT(solver().isTrue(~b));
		solver().propagate();
		solver().assume(~lits[0].first) && solver().propagate();
		CPPUNIT_ASSERT(solver().isTrue(lits[1].first));
	}

	void testTrivialBackpropFalse() {
		LitVec lits = makeLits();
		solver().force(~body, 0);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 1));
		CPPUNIT_ASSERT(ctx.numConstraints() == 0);
		CPPUNIT_ASSERT(solver().isFalse(lits[1]));
		CPPUNIT_ASSERT(solver().isFalse(lits[2]));
		CPPUNIT_ASSERT(solver().isFalse(lits[3]));
		CPPUNIT_ASSERT(solver().isFalse(lits[4]));
	}

	void testTrivialBackpropFalseWeight() {
		WeightLitVec lits = makeWeightLits();
		solver().force(~body, 0);
		CPPUNIT_ASSERT_EQUAL(true,  WeightConstraint::newWeightConstraint(ctx, body, lits, 2));
		CPPUNIT_ASSERT(ctx.numConstraints() == 1);
		CPPUNIT_ASSERT(solver().isFalse(a));
		CPPUNIT_ASSERT(solver().isFalse(~b));
	}

	void testForwardTrue() {
		LitVec assume, expected;
		assume.push_back(a);
		assume.push_back(~c);
		expected.push_back(body);
		propCard(assume, expected);
	}

	void testForwardFalse() {
		LitVec assume;
		assume.push_back(~a);
		assume.push_back(c);
		assume.push_back(~d);
		
		LitVec expect;
		expect.push_back(~body);

		propCard(assume, expect);
	}

	void testBackwardTrue() {
		LitVec assume, expect;
		assume.push_back(body);
		assume.push_back(c);
		assume.push_back(~d);
		expect.push_back(a);
		expect.push_back(~b);
		propCard(assume, expect);
		
	}

	void testBackwardFalse() {
		LitVec assume, expect;
		assume.push_back(~body);
		assume.push_back(d);
		expect.push_back( ~a );
		expect.push_back( b );
		expect.push_back( c );
		propCard(assume, expect);
	}

	void testForwardTrueConflict() {
		LitVec assume;
		assume.push_back(a);
		assume.push_back(~c);
		propConflictTest(assume, ~body);
	}

	void testForwardFalseConflict() {
		LitVec assume;
		assume.push_back(~a);
		assume.push_back(c);
		assume.push_back(~d);
		propConflictTest(assume, body);
		
	}

	void testBackwardTrueConflict() {
		LitVec assume;
		assume.push_back(body);
		assume.push_back(c);
		assume.push_back(~d);
		propConflictTest(assume, b);
	}

	void testBackwardFalseConflict() {
		LitVec assume, expect;
		assume.push_back(~body);
		assume.push_back(d);
		propConflictTest(assume, ~b);
	}

	void testReasonBug() {
		LitVec lits = makeLits();
		lits.push_back(~e);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 3));
		LitVec assume, reason;
		assume.push_back(a);
		assume.push_back(~b);
		assume.push_back(~d);
		assume.push_back(e);
		for (uint32 i = 0; i < assume.size(); ++i) {
			CPPUNIT_ASSERT_EQUAL(true, solver().assume(assume[i]));
			CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		}
		CPPUNIT_ASSERT(assume.size() == solver().numAssignedVars());
		
		// B -> ~c because of: ~d, e, B
		CPPUNIT_ASSERT_EQUAL(true, solver().assume(body));
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(~c));
		//CPPUNIT_ASSERT(con == solver().reason(c.var()).constraint());
		solver().reason(~c, reason);
		CPPUNIT_ASSERT_EQUAL(LitVec::size_type(3), reason.size());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), ~d) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), e) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), body) != reason.end());
		solver().undoUntil(solver().decisionLevel()-1);
		reason.clear();

		// ~B -> c because of: a, ~b, ~B
		CPPUNIT_ASSERT_EQUAL(true, solver().assume(~body));
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(c));
		// CPPUNIT_ASSERT(con == solver().reason(c.var()).constraint());
		solver().reason(c, reason);
		CPPUNIT_ASSERT_EQUAL(LitVec::size_type(3), reason.size());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), a) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), ~b) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), ~body) != reason.end());
		solver().undoUntil(solver().decisionLevel()-1);
		reason.clear();

		// ~c -> B because of: a, ~b, ~c
		CPPUNIT_ASSERT_EQUAL(true, solver().assume(~c));
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(body));
		//CPPUNIT_ASSERT(con == solver().reason(body.var()).constraint());
		solver().reason(body, reason);
		CPPUNIT_ASSERT_EQUAL(LitVec::size_type(3), reason.size());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), a) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), ~b) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), ~c) != reason.end());
		solver().undoUntil(solver().decisionLevel()-1);

		// c -> ~B because of: ~d, e, c
		CPPUNIT_ASSERT_EQUAL(true, solver().assume(c));
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(~body));
		//CPPUNIT_ASSERT(con == solver().reason(body.var()).constraint());
		solver().reason(~body, reason);
		CPPUNIT_ASSERT_EQUAL(LitVec::size_type(3), reason.size());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), ~d) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), e) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), c) != reason.end());
		solver().undoUntil(solver().decisionLevel()-1);
		reason.clear();
	}

	void testWeightReasonAfterBackprop() {
		WeightLitVec lits = makeWeightLits();
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		solver().assume(~body) && solver().propagate();
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(~a));
		solver().assume(d) && solver().propagate();
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(b));
		LitVec r;
		solver().reason(~a, r);
		CPPUNIT_ASSERT(r.size() == 1 && r[0] == ~body);
		solver().reason(b, r);
		CPPUNIT_ASSERT(r.size() == 2 && r[0] == ~body && r[1] == d);
		solver().undoUntil(solver().decisionLevel()-1);
		solver().reason(~a, r);
		CPPUNIT_ASSERT(r.size() == 1 && r[0] == ~body);
		solver().undoUntil(solver().decisionLevel()-1);
	}
	
	void testOrderBug() {
		LitVec lits;
		lits.push_back(body);
		lits.push_back(a);
		lits.push_back(b);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 1));
		solver().assume(e) && solver().propagate();
		
		solver().force(~a, 0);
		solver().force(body, 0);
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(b));
		LitVec reason;
		solver().reason(b, reason);
		CPPUNIT_ASSERT(LitVec::size_type(2) == reason.size());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), body) != reason.end());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), ~a) != reason.end());

	}

	void testBackwardAfterForward() {
		LitVec lits;
		lits.push_back(body);
		lits.push_back(a);
		lits.push_back(b);
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 1));
		solver().assume(a);
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(body));
		LitVec reason;
		solver().reason(body, reason);
		CPPUNIT_ASSERT(LitVec::size_type(1) == reason.size());
		CPPUNIT_ASSERT(std::find(reason.begin(), reason.end(), a) != reason.end());

		solver().assume(~b);
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(body));
	}

	void testSimplify() {
		LitVec lits = makeLits();
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 2));
		CPPUNIT_ASSERT_EQUAL(true, solver().simplify());
		CPPUNIT_ASSERT_EQUAL(1u, ctx.numConstraints());
		solver().force(a, 0);
		solver().simplify();
		CPPUNIT_ASSERT_EQUAL(1u, ctx.numConstraints());
		solver().force(~c, 0);
		solver().simplify();
		CPPUNIT_ASSERT_EQUAL(0u, ctx.numConstraints());
	}

	void testSimplifyCardinality() {
		LitVec lits = makeLits();
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 2));
		ctx.addUnary(~a);
		ctx.addUnary(~d);
		ctx.addUnary(~b);
		solver().simplify();
		solver().assume(~c);
		solver().propagate();
		CPPUNIT_ASSERT(solver().isTrue(body));
		LitVec out;
		solver().reason(body, out);
		// a, d, and ~b were removed from constraint
		CPPUNIT_ASSERT(out.size() == 1 && out[0] == ~c);
	}
	void testSimplifyWeight() {
		WeightLitVec lits = makeWeightLits();
		lits[2].second    = 2;
		lits.push_back(WeightLiteral(e, 1));
		lits.push_back(WeightLiteral(f, 1));
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 2));
		ctx.addUnary(body);
		ctx.addUnary(~a);
		ctx.addUnary(~d);
		ctx.addUnary(~e);
		ctx.addUnary(~f);
		solver().simplify();
		CPPUNIT_ASSERT(solver().value(b.var()) == value_free);
		solver().assume(c);
		solver().propagate();
		CPPUNIT_ASSERT(solver().isTrue(~b));
		LitVec out;
		solver().reason(~b, out);
		CPPUNIT_ASSERT(out.size() == 1 && out[0] == c);
	}

	void testAssertWeightTriviallySat() {
		WeightLitVec lits;
		lits.push_back(WeightLiteral(a, 2));
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 0));
		CPPUNIT_ASSERT(solver().isTrue(body));
	}
	void testAssertWeightTriviallyUnSat() {
		WeightLitVec lits;
		lits.push_back(WeightLiteral(a, 2));
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		CPPUNIT_ASSERT(solver().isFalse(body));
	}

	void testAssertWeightNotSoTriviallySat() {
		WeightLitVec lits = makeWeightLits();
		solver().force(lits[1].first, 0);
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body,lits, 2));
		CPPUNIT_ASSERT(solver().isTrue(body));
	}

	void testAssertWeightNotSoTriviallyUnSat() {
		WeightLitVec lits = makeWeightLits();
		solver().force(~lits[0].first, 0);
		solver().force(~lits[2].first, 0);
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 4));
		CPPUNIT_ASSERT(solver().isTrue(~body));
	}

	void testWeightForwardTrue() {
		LitVec assume, expect;
		assume.push_back(a);
		expect.push_back(body);
		propWeight(assume, expect);   

		assume.clear();
		assume.push_back(~b);
		assume.push_back(~c);
		propImpl(assume, expect);

		assume.clear();
		assume.push_back(~b);
		assume.push_back(d);
		propImpl(assume, expect);
	}

	void testWeightForwardFalse() {
		LitVec assume, expect;
		assume.push_back(~a);
		assume.push_back(b);
		expect.push_back(~body);
		propWeight(assume, expect);   
	}

	void testWeightBackwardTrue() {
		WeightLitVec lits = makeWeightLits();
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		solver().assume(~a);
		solver().force(body, 0);
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(~b));
		CPPUNIT_ASSERT_EQUAL(value_free, solver().value(c.var()));
		LitVec r;
		solver().reason(~b, r);
		CPPUNIT_ASSERT(r.size() == 2);
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), ~a) != r.end());
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), body) != r.end());
		
		solver().assume(~d);
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(~c));

		solver().reason(~c, r);
		CPPUNIT_ASSERT(r.size() == 3);
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), ~a) != r.end());
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), body) != r.end());
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), ~d) != r.end());

		solver().undoUntil(solver().decisionLevel()-1);
		solver().reason(~b, r);
		CPPUNIT_ASSERT(r.size() == 2);
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), ~a) != r.end());
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), body) != r.end());
	}

	void testWeightBackwardFalse() {
		WeightLitVec lits = makeWeightLits();
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		solver().assume(~body);
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(~a));
		LitVec r;
		solver().reason(~a, r);
		CPPUNIT_ASSERT(r.size() == 1);
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), ~body) != r.end());
		
		solver().force(~b, 0);
		CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(c));
		CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(~d));

		LitVec r2;
		solver().reason(c, r);
		solver().reason(~d, r2);
		CPPUNIT_ASSERT(r == r2);
		CPPUNIT_ASSERT(r.size() == 2);
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), ~body) != r.end());
		CPPUNIT_ASSERT(std::find(r.begin(), r.end(), ~b) != r.end());
	}

	void testWeightConflict() {
		WeightLitVec lits = makeWeightLits();
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		LitVec assume;
		assume.push_back(body);
		assume.push_back(~a);
		assume.push_back(b);
		std::sort(assume.begin(), assume.end());
		do {
			CPPUNIT_ASSERT_EQUAL(true, solver().assume(assume[0]));
			for (uint32 i = 1; i < assume.size(); ++i) {  
				CPPUNIT_ASSERT_EQUAL(true, solver().force(assume[i],0));
			}
			CPPUNIT_ASSERT_EQUAL(false, solver().propagate());
			solver().undoUntil(0);
		} while (std::next_permutation(assume.begin(), assume.end()));
	}

	void testCloneWeight() {
		ctx.setSolvers(2);
		WeightLitVec lits = makeWeightLits();
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		solver().force(~a, 0);
		ctx.endInit();
		Solver solver2;
		ctx.attach(solver2);
		CPPUNIT_ASSERT(solver2.numConstraints() == 1);
		
		CPPUNIT_ASSERT(solver2.numWatches(a) == 0 && solver2.numWatches(~a) == 0);
		solver2.assume(body);
		solver2.propagate();
		solver2.assume(~d);
		solver2.propagate();
		CPPUNIT_ASSERT(solver2.isTrue(~b));
		CPPUNIT_ASSERT(solver2.isTrue(~c));
		LitVec out;
		solver2.reason(~b, out);
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), body) != out.end());
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), ~a) != out.end());
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), ~d) == out.end());
		out.clear();
		solver2.reason(~c, out);
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), ~d) != out.end());
	}

	void testCloneWeightShared() {
		ctx.setSolvers(2);
		ctx.physicalSharing(SharedContext::share_problem);
		WeightLitVec lits = makeWeightLits();
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		solver().force(~a, 0);
		ctx.endInit();
		Solver solver2;
		ctx.attach(solver2);
		CPPUNIT_ASSERT(solver2.numConstraints() == 1);
		
		CPPUNIT_ASSERT(solver2.numWatches(a) == 0 && solver2.numWatches(~a) == 0);
		solver2.assume(body);
		solver2.propagate();
		solver2.assume(~d);
		solver2.propagate();
		CPPUNIT_ASSERT(solver2.isTrue(~b));
		CPPUNIT_ASSERT(solver2.isTrue(~c));
		LitVec out;
		solver2.reason(~b, out);
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), body) != out.end());
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), ~a) != out.end());
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), ~d) == out.end());
		out.clear();
		solver2.reason(~c, out);
		CPPUNIT_ASSERT(std::find(out.begin(), out.end(), ~d) != out.end());
	}
private:
	static bool newCardinalityConstraint(SharedContext& ctx, const LitVec& lits, int bound) {
		CPPUNIT_ASSERT(lits.size() >  1);
		WeightLitVec wlits;
		for (LitVec::size_type i = 1; i < lits.size(); ++i) {
			wlits.push_back(WeightLiteral(lits[i], 1));
		}
		return WeightConstraint::newWeightConstraint(ctx, lits[0], wlits, bound);
	}
	void propCard(LitVec& assumptions, const LitVec& expected) {
		LitVec lits = makeLits();
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 2));
		propImpl(assumptions, expected);
	}
	void propWeight(LitVec& assume, LitVec& expect) {
		WeightLitVec lits = makeWeightLits();
		CPPUNIT_ASSERT_EQUAL(true, WeightConstraint::newWeightConstraint(ctx, body, lits, 3));
		propImpl(assume, expect);
	}
	void propImpl(LitVec& assumptions, const LitVec& expected) {
		std::sort(assumptions.begin(), assumptions.end());
		do {
			for (uint32 i = 0; i < assumptions.size(); ++i) {
				CPPUNIT_ASSERT_EQUAL(true, solver().assume(assumptions[i]));
				CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
			}
			for (uint32 i = 0; i < expected.size(); ++i) {
				CPPUNIT_ASSERT_EQUAL(true, solver().isTrue(expected[i]));
				LitVec reason;
				solver().reason(expected[i], reason);
				CPPUNIT_ASSERT_EQUAL(assumptions.size(), reason.size());
				for (uint32 j = 0; j < assumptions.size(); ++j) {
					CPPUNIT_ASSERT(find(reason.begin(), reason.end(), assumptions[j]) != reason.end());
				}
				
			}
			solver().undoUntil(0);
		} while (std::next_permutation(assumptions.begin(), assumptions.end()));
	}
	void propConflictTest(LitVec& assumptions, Literal cflLit) {
		LitVec lits = makeLits();
		CPPUNIT_ASSERT_EQUAL(true, newCardinalityConstraint(ctx, lits, 2));
		do {
			for (uint32 i = 0; i < assumptions.size()-1; ++i) {
				CPPUNIT_ASSERT_EQUAL(true, solver().assume(assumptions[i]));
				CPPUNIT_ASSERT_EQUAL(true, solver().propagate());
			}
			CPPUNIT_ASSERT_EQUAL(true, solver().assume(assumptions.back()));
			CPPUNIT_ASSERT_EQUAL(true, solver().force(cflLit, 0));
			CPPUNIT_ASSERT_EQUAL(false, solver().propagate());
			LitVec cfl = solver().conflict();
			CPPUNIT_ASSERT_EQUAL(assumptions.size() + 1, cfl.size());
			for (uint32 i = 0; i < assumptions.size(); ++i) {
				CPPUNIT_ASSERT(std::find(cfl.begin(), cfl.end(), assumptions[i]) != cfl.end());
			}
			CPPUNIT_ASSERT(std::find(cfl.begin(), cfl.end(), cflLit) != cfl.end());
			solver().undoUntil(0);
		}while (std::next_permutation(assumptions.begin(), assumptions.end()));
	} 
	
	LitVec makeLits() {
		LitVec res;
		res.push_back(body);
		res.push_back(a);
		res.push_back(~b);
		res.push_back(~c);
		res.push_back(d);
		return res;
	}
	WeightLitVec makeWeightLits() {
		WeightLitVec res;
		res.push_back(WeightLiteral(a, 4));
		res.push_back(WeightLiteral(~b, 2));
		res.push_back(WeightLiteral(~c, 1));
		res.push_back(WeightLiteral(d, 1));
		return res;
	}
	SharedContext ctx;
	Literal body;
	Literal a, b, c, d, e, f;
};
CPPUNIT_TEST_SUITE_REGISTRATION(CardTest);
} } 
