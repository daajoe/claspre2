// 
// Copyright (c) 2006-2012, Benjamin Kaufmann
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
#include <clasp/solver.h>
#include <clasp/clause.h>

#ifdef _MSC_VER
#pragma warning (disable : 4996) // 'std::copy': Function call with parameters that may be unsafe
#endif

namespace Clasp { 

DecisionHeuristic::~DecisionHeuristic() {}
/////////////////////////////////////////////////////////////////////////////////////////
// SelectFirst selection strategy
/////////////////////////////////////////////////////////////////////////////////////////
// selects the first free literal
Literal SelectFirst::doSelect(Solver& s) {
	for (Var i = 1; i <= s.numVars(); ++i) {
		if (s.value(i) == value_free) {
			return selectLiteral(s, i, 0);
		}
	}
	assert(!"SelectFirst::doSelect() - precondition violated!\n");
	return Literal();
}
/////////////////////////////////////////////////////////////////////////////////////////
// Post propagator list
/////////////////////////////////////////////////////////////////////////////////////////
Solver::PPList::PPList() : head(0), look(0) { }
Solver::PPList::~PPList() {
	enableFunky();
	for (PostPropagator* r = head; r;) {
		PostPropagator* t = r;
		r = r->next;
		t->destroy();
	}
}
PostPropagator* Solver::PPList::disableFunky() {
	// already disabled
	if (look >= nil()) { return look != nil() ? look : 0; }
	// move post propagators with priority >= lookahead to ext
	PostPropagator* ext = 0;
	if (head && head->priority() < PostPropagator::priority_lookahead) {
		for (PostPropagator* r = head; r->next; r = r->next) {
			if (r->next->priority() >= PostPropagator::priority_lookahead) {
				ext     = r->next;
				r->next = 0;
				break;
			}
		}
	}
	else {
		ext  = head;
		head = 0;
	}
	// mark as disabled even if there is currently no funky post propagator
	look = ext ? ext : nil();
	return ext;
}
PostPropagator* Solver::PPList::enableFunky() {
	// already enabled or empty list
	if (look <= nil()) { return (look = 0); }
	// move post propagators with priority >= lookahead back to active list
	PostPropagator* ext = look;
	look                = 0; // mark as enabled
	PostPropagator* r  = head;
	if (!r)         { return (head = ext); }
	while (r->next) { r= r->next; }
	return (r->next = ext);	
}

void Solver::PPList::add(PostPropagator* p) {
	assert(p && p->next == 0);
	uint32 prio = p->priority();
	PostPropagator** listHead = (look == 0 || prio < PostPropagator::priority_lookahead)
		? &head
		: &look;
	if (*listHead <= nil() || prio < (*listHead)->priority()) {
		p->next  = *listHead > nil() ? *listHead : 0;
		*listHead= p;
		return;
	}
	for (PostPropagator* r = *listHead;; r = r->next) {
		if (!r->next || prio < r->next->priority()) {
			p->next = r->next;
			r->next = p;
			return;
		}
	}
}

void Solver::PPList::remove(PostPropagator* p) {
	assert(p);
	PostPropagator** listHead = (look == 0 || p->priority() < PostPropagator::priority_lookahead)
		? &head
		: &look;
	if (*listHead <= nil()) { return; }
	if (*listHead == p)     {
		*listHead = (*listHead)->next;
		p->next   = 0;
		return;
	}
	for (PostPropagator* r = *listHead; r->next; r = r->next) {
		if (r->next == p) {
			r->next = r->next->next;
			p->next = 0;
			return;
		}
	}
}

bool Solver::PPList::propagate(Solver& s, PostPropagator* x) const {
	for (PostPropagator* p = head, *t; p != x;) {
		// just in case t removes itself from the list
		// during propagateFixpoint
		t = p;
		p = p->next;
		if (!t->propagateFixpoint(s)) { return false; }
	}
	return true;
}
void Solver::PPList::simplify(Solver& s, bool shuf) {
	for (PostPropagator* r = head; r;) {
		PostPropagator* t = r;
		r = r->next;
		if (t->simplify(s, shuf)) {
			remove(t);
		}
	}
}
bool Solver::PPList::init(Solver& s, PostPropagator* list) const {
	for (PostPropagator* r = list; r; r = r->next) {
		if (!r->init(s)) return false;
	}
	return true;
}
bool PostPropagator::propagateFixpoint(Solver& s) {
	bool ok = propagate(s);
	while (ok && s.queueSize() > 0) {
		ok = s.propagateUntil(this) && propagate(s);
	}
	return ok;
}
void Solver::PPList::reset() const      { for (PostPropagator* r = head; r; r = r->next) { r->reset(); } }
bool Solver::PPList::isModel(Solver& s) const {
	if (s.hasConflict()) { return false; }
	for (PostPropagator* r = head; r; r = r->next) {
		if (!r->isModel(s)) return false;
	}
	return true;
}
/////////////////////////////////////////////////////////////////////////////////////////
// SolverStrategies
/////////////////////////////////////////////////////////////////////////////////////////
SolverStrategies::SolverStrategies() {
	struct X { RNG y; uint32 z[3]; };
	static_assert(sizeof(SolverStrategies) == sizeof(X), "Unsupported Padding");
	std::memset(&compress, 0, 3*sizeof(uint32));
	ccMinAntes = all_antes;
	search     = use_learning;
}
SolverStrategies::HeuFactory SolverStrategies::heuFactory_s = 0;
/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Construction/Destruction/Setup
////////////////////////////////////////////////////////////////////////////////////////
Solver::Solver() 
	: strategy_()
	, ccMin_(0)
	, smallAlloc_(new SmallClauseAlloc)
	, shared_(0)
	, undoHead_(0)
	, enum_(0)
	, memLimit_(0)
	, id_(0)
	, units_(0)
	, lastSimplify_(0)
	, rootLevel_(0)
	, btLevel_(0)
	, lbdTime_(0)
	, shuffle_(false) {
	Var sentVar = assign_.addVar();
	assign_.setValue(sentVar, value_true);
	markSeen(sentVar);
	setMemLimit( 16384 ); // 16 GB
}

Solver::~Solver() {
	freeMem();
}

void Solver::freeMem() {
	std::for_each( constraints_.begin(), constraints_.end(), DestroyObject());
	std::for_each( learnts_.begin(), learnts_.end(), DestroyObject() );
	constraints_.clear();
	learnts_.clear();
	setEnumerationConstraint(0);
	PodVector<WatchList>::destruct(watches_);
	// free undo lists
	// first those still in use
	for (DecisionLevels::size_type i = 0; i != levels_.size(); ++i) {
		delete levels_[i].undo;
	}
	// then those in the free list
	for (ConstraintDB* x = undoHead_; x; ) {
		ConstraintDB* t = x;
		x = (ConstraintDB*)x->front();
		delete t;
	}
	delete smallAlloc_;
	delete ccMin_;
	smallAlloc_   = 0;
	ccMin_        = 0;
}

SatPreprocessor* Solver::satPrepro() const { return shared_->satPrepro.get(); }

void Solver::reset() {
	// hopefully, no one derived from this class...
	this->~Solver();
	new (this) Solver();
}

void Solver::setMemLimit(uint32 learntLimitInMB) {
	memLimit_   = learntLimitInMB;
	memLimit_ <<= 20;
}

void Solver::startInit(SharedContext* ctx, uint32 numConsGuess) {
	shared_ = ctx;
	assert(numVars() <= shared_->numVars());
	assign_.resize(shared_->numVars() + 1);
	watches_.resize(assign_.numVars()<<1);
	// pre-allocate some memory
	assign_.trail.reserve(numVars());
	constraints_.reserve(numConsGuess/2);
	levels_.reserve(25);
	if (smallAlloc_ == 0) { smallAlloc_ = new SmallClauseAlloc(); }
	if (undoHead_ == 0)   {
		for (uint32 i = 0; i != 25; ++i) { 
			undoFree(new ConstraintDB(10)); 
		}
	}
	if (strategy_.compress == 0) { strategy_.compress = UINT32_MAX; }
	if (heuristic_.get() == 0)   {
		if (SolverStrategies::heuFactory_s == 0) { strategy_.heuId = 0; }
		if (strategy_.heuId == 0) { heuristic_.reset(new SelectFirst()); }
		else                      { heuristic_.reset(strategy_.heuFactory_s(strategy_)); }
	}
	heuristic_->startInit(*this);
	// disable any funky post propagators during setup
	post_.disableFunky();
}

bool Solver::endInit() {
	if (!propagate()) { return false; }
	heuristic_->endInit(*this);
	if (strategy_.signFix) {
		for (Var v = 1; v <= numVars(); ++v) {
			if (savedValue(v) != value_free) { initPrefValue(v, valSign(savedValue(v))); }
			else                             { initPrefValue(v, 1 + defaultLiteral(v).sign()); }
		}
	}
	// enable funky post propagators and
	// force initial propagation
	PostPropagator* ext = post_.enableFunky();
	if (ext && (!post_.init(*this, ext) || !propagate())) {
		return false;
	}
	return simplify();
}

void Solver::add(Constraint* c) {
	assert(!shared_->frozen() && this == shared_->master());
	constraints_.push_back(c);
}

bool Solver::addUnary(Literal p, ConstraintType t) {
	assert(t != Constraint_t::static_constraint || (!shared_->frozen() && this == shared_->master()));
	if (t != Constraint_t::static_constraint) {
		bool subsumed = (isTrue(p) && level(p.var()) == 0) || (shared_->isShared() && shared_->allowImplicit(t) && !shared_->addBinary(p, negLit(0), true));
		if (!subsumed) {
			stats.addLearnt(1, t);
			shared_->distribute(*this, &p, 1, ClauseInfo(t).setLbd(1));
		}
	}
	return force(p, 0, Antecedent(posLit(0)));
}

void Solver::addShort(const Literal* c, uint32 sz, const ClauseInfo& x) {
	assert(sz == 2 || sz == 3);
	bool learnt  = x.learnt();
	if (!learnt || shared_->allowImplicit(x.type())) {
		learnt = (sz == 3 ? shared_->addTernary(c[0], c[1], c[2], learnt) : shared_->addBinary(c[0], c[1], learnt)) && learnt;
		if (!learnt) { return; }
	}
	else { learnts_.push_back(Clause::newClause(*this, c, sz, x)); }
	stats.addLearnt(sz, x.type());
	shared_->distribute(*this, c, sz, x);
}

uint32 Solver::numConstraints() const {
	return static_cast<uint32>(constraints_.size())
		+ (shared_ ? shared_->numBinary()+shared_->numTernary() : 0);
}

bool Solver::popRootLevel(uint32 i)  {
	clearStopConflict();
	if (i > rootLevel_) i = rootLevel_;
	rootLevel_ = btLevel_ = rootLevel_ - i;
	impliedLits_.front    = 0;
	// go back to new root level and re-assert still implied literals
	undoUntil(rootLevel_, true);
	if (!isTrue(sharedContext()->tagLiteral())) {
		removeConditional();
	}
	return !hasConflict();
}

bool Solver::clearAssumptions()  {
	return popRootLevel(rootLevel())
		&& simplify();
}

void Solver::clearStopConflict() {
	if (hasConflict() && hasStopConflict()) {
		rootLevel_ = conflict_[1].asUint();
		conflict_.clear();
	}
}

bool Solver::hasStopConflict() const { return conflict_[0] == negLit(0); }
void Solver::setStopConflict() {
	if (!hasConflict()) {
		// we use the nogood {FALSE} to represent the unrecoverable conflict -
		// note that {FALSE} can otherwise never be a violated nogood because
		// TRUE is always true in every solver
		conflict_.push_back(negLit(0));
		// remember the current root-level
		conflict_.push_back(Literal::fromRep(rootLevel_));
	}
	// artificially increase root level -
	// this way, the solver is prevented from resolving the conflict
	pushRootLevel(decisionLevel());
}

void Solver::updateGuidingPath(LitVec& gpOut, LitVec::size_type& start, uint32& implied) {
	if (start < shared_->topLevelSize()) {
		start = shared_->topLevelSize();
	}
	const LitVec& tr = assign_.trail;
	LitVec::size_type end = rootLevel_ == decisionLevel() ? tr.size() : levels_[rootLevel_].trailPos;
	for (LitVec::size_type i = start; i < end; ++i) {
		if (reason(tr[i]).isNull()) {
			gpOut.push_back(tr[i]);
		}
	}
	start = end;
	uint32 implLits = 0;
	for (ImpliedList::iterator it = impliedLits_.begin(); it != impliedLits_.end(); ++it) {
		if (it->level <= rootLevel_ && ++implLits > implied) {
			gpOut.push_back(it->lit);
		}
	}
	implied = implLits;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Watch management
////////////////////////////////////////////////////////////////////////////////////////
uint32 Solver::numWatches(Literal p) const {
	assert( validVar(p.var()) );
	if (!validWatch(p)) return 0;
	return static_cast<uint32>(watches_[p.index()].size()) 
		+ shared_->numShortImplications(p);
}
	
bool Solver::hasWatch(Literal p, Constraint* c) const {
	if (!validWatch(p)) return false;
	const WatchList& pList = watches_[p.index()];
	return std::find_if(pList.right_begin(), pList.right_end(), GenericWatch::EqConstraint(c)) != pList.right_end();
}

bool Solver::hasWatch(Literal p, ClauseHead* h) const {
	if (!validWatch(p)) return false;
	const WatchList& pList = watches_[p.index()];
	return std::find_if(pList.left_begin(), pList.left_end(), ClauseWatch::EqHead(h)) != pList.left_end();
}

GenericWatch* Solver::getWatch(Literal p, Constraint* c) const {
	if (!validWatch(p)) return 0;
	const WatchList& pList = watches_[p.index()];
	WatchList::const_right_iterator it = std::find_if(pList.right_begin(), pList.right_end(), GenericWatch::EqConstraint(c));
	return it != pList.right_end()
		? &const_cast<GenericWatch&>(*it)
		: 0;
}

void Solver::removeWatch(const Literal& p, Constraint* c) {
	assert(validWatch(p));
	WatchList& pList = watches_[p.index()];
	pList.erase_right(std::find_if(pList.right_begin(), pList.right_end(), GenericWatch::EqConstraint(c)));
}

void Solver::removeWatch(const Literal& p, ClauseHead* h) {
	assert(validWatch(p));
	WatchList& pList = watches_[p.index()];
	pList.erase_left(std::find_if(pList.left_begin(), pList.left_end(), ClauseWatch::EqHead(h)));
}

bool Solver::removeUndoWatch(uint32 dl, Constraint* c) {
	assert(dl != 0 && dl <= decisionLevel() );
	if (levels_[dl-1].undo) {
		ConstraintDB& uList = *levels_[dl-1].undo;
		ConstraintDB::iterator it = std::find(uList.begin(), uList.end(), c);
		if (it != uList.end()) {
			*it = uList.back();
			uList.pop_back();
			return true;
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Basic DPLL-functions
////////////////////////////////////////////////////////////////////////////////////////

// removes all satisfied binary and ternary clauses as well
// as all constraints for which Constraint::simplify returned true.
bool Solver::simplify() {
	if (decisionLevel() != 0) return true;
	if (hasConflict())        return false;
	if (lastSimplify_ != assign_.trail.size()) {
		LitVec::size_type old = lastSimplify_;
		if (!simplifySAT()) { return false; }
		assert(lastSimplify_ == assign_.trail.size());
		heuristic_->simplify(*this, old);
	}
	if (shuffle_) { simplifySAT(); }
	return true;
}

void Solver::removeConditional() { 
	Literal p = ~sharedContext()->tagLiteral();
	if (!isSentinel(p)) {
		ConstraintDB::size_type i, j, end = learnts_.size();
		for (i = j = 0; i != end; ++i) {
			ClauseHead* c = static_cast<LearntConstraint*>(learnts_[i])->clause();
			if (!c || !c->tagged()) {
				learnts_[j++] = c;
			}
			else {
				assert((decisionLevel() == rootLevel() || !c->locked(*this)) && "Solver::removeConditional(): must not remove locked constraint!");
				c->destroy(this, true);
			}
		}
		learnts_.erase(learnts_.begin()+j, learnts_.end());
	}
}

void Solver::strengthenConditional() { 
	Literal p = ~sharedContext()->tagLiteral();
	if (!isSentinel(p)) {
		ConstraintDB::size_type i, j, end = learnts_.size();
		for (i = j = 0; i != end; ++i) {
			ClauseHead* c = static_cast<LearntConstraint*>(learnts_[i])->clause();
			if (!c || !c->tagged() || !c->strengthen(*this, p, true).second) {
				learnts_[j++] = c;
			}
			else {
				assert((decisionLevel() == rootLevel() || !c->locked(*this)) && "Solver::strengthenConditional(): must not remove locked constraint!");
				c->destroy(this, false);
			}
		}
		learnts_.erase(learnts_.begin()+j, learnts_.end());
	}
}

bool Solver::simplifySAT() {
	if (queueSize() > 0 && !propagate()) {
		return false;
	}
	assert(assign_.qEmpty());
	assign_.front = lastSimplify_;
	while (!assign_.qEmpty()) {
		simplifyShort(assign_.qPop()); // remove satisfied binary- and ternary clauses
	}
	lastSimplify_ = assign_.front;
	if (shuffle_) {
		std::random_shuffle(constraints_.begin(), constraints_.end(), strategies().rng);
		std::random_shuffle(learnts_.begin(), learnts_.end(), strategies().rng);
	}
	simplifyDB(constraints_);
	simplifyDB(learnts_);
	post_.simplify(*this, shuffle_);
	if (enum_ && enum_->simplify(*this, shuffle_)) {
		enum_->destroy(this, false);
		enum_ = 0;
	}
	shuffle_ = false;
	return true;
}

void Solver::simplifyDB(ConstraintDB& db) {
	ConstraintDB::size_type i, j, end = db.size();
	for (i = j = 0; i != end; ++i) {
		Constraint* c = db[i];
		if (c->simplify(*this, shuffle_)) { c->destroy(this, false); }
		else                              { db[j++] = c;  }
	}
	db.erase(db.begin()+j, db.end());
}

// removes all binary clauses containing p - those are now SAT
// binary clauses containing ~p are unit and therefore likewise SAT. Those
// are removed when their second literal is processed.
// Note: Binary clauses containing p are those that watch ~p.
//
// Simplifies ternary clauses.
// Ternary clauses containing p are SAT and therefore removed.
// Ternary clauses containing ~p are now either binary or SAT. Those that
// are SAT are removed when the satisfied literal is processed. 
// All conditional binary-clauses are replaced with a real binary clause.
// Note: Ternary clauses containing p watch ~p. Those containing ~p watch p.
// Note: Those clauses are now either binary or satisfied.
void Solver::simplifyShort(Literal p) {
	releaseVec(watches_[p.index()]);
	releaseVec(watches_[(~p).index()]);
	if (!shared_->isShared()) {
		shared_->shortImplications().removeTrue(*this, p);
	}
}

void Solver::setConflict(Literal p, const Antecedent& a, uint32 data) {
	++stats.conflicts;
	conflict_.push_back(~p);
	if (strategy_.search != SolverStrategies::no_learning && !a.isNull()) {
		if (data == UINT32_MAX) {
			a.reason(*this, p, conflict_);
		}
		else {
			// temporarily replace old data with new data
			uint32 saved = assign_.data(p.var());
			assign_.setData(p.var(), data);
			// extract conflict using new data
			a.reason(*this, p, conflict_);
			// restore old data
			assign_.setData(p.var(), saved);
		}
	}
}

bool Solver::assume(const Literal& p) {
	if (value(p.var()) == value_free) {
		assert(decisionLevel() != assign_.maxLevel());
		++stats.choices;
		levels_.push_back(DLevel(numAssignedVars(), 0));
		return assign_.assign(p, decisionLevel(), Antecedent());
	}
	return isTrue(p);
}

bool Solver::propagate() {
	if (unitPropagate() && post_.propagate(*this, 0)) {
		assert(queueSize() == 0);
		return true;
	}
	assign_.qReset();
	post_.reset();
	return false;
}

Constraint::PropResult ClauseHead::propagate(Solver& s, Literal p, uint32&) {
	Literal* head = head_;
	uint32 wLit   = (head[1] == ~p); // pos of false watched literal
	if (s.isTrue(head[1-wLit])) {
		return Constraint::PropResult(true, true);
	}
	else if (!s.isFalse(head[2])) {
		assert(!isSentinel(head[2]) && "Invalid ClauseHead!");
		head[wLit] = head[2];
		head[2]    = ~p;
		s.addWatch(~head[wLit], ClauseWatch(this));
		return Constraint::PropResult(true, false);
	}
	else if (updateWatch(s, wLit)) {
		assert(!s.isFalse(head_[wLit]));
		s.addWatch(~head[wLit], ClauseWatch(this));
		return Constraint::PropResult(true, false);
	}	
	return PropResult(s.force(head_[1^wLit], this), true);
}

bool Solver::propagateUnits() {
	assert(decisionLevel() == 0 && queueSize() == 0);
	assign_.front = units_;
	units_        = (uint32)assign_.trail.size();
	while (!assign_.qEmpty()) { markSeen(assign_.qPop().var()); } 
	return true;
}

bool Solver::unitPropagate() {
	assert(!hasConflict());
	Literal p, q, r;
	uint32 idx, ignore, DL = decisionLevel();
	const ShortImplicationsGraph& btig = shared_->shortImplications();
	while ( !assign_.qEmpty() ) {
		p             = assign_.qPop();
		idx           = p.index();
		WatchList& wl = watches_[idx];
		// first: short clause BCP
		if (!btig.propagate(*this, p)) {
			return false;
		}
		// second: clause BCP
		if (wl.left_size() != 0) {
			WatchList::left_iterator it, end, j = wl.left_begin(); 
			Constraint::PropResult res;
			for (it = wl.left_begin(), end = wl.left_end();  it != end;  ) {
				ClauseWatch& w = *it++;
				res = w.head->ClauseHead::propagate(*this, p, ignore);
				if (res.keepWatch) {
					*j++ = w;
				}
				if (!res.ok) {
					wl.shrink_left(std::copy(it, end, j));
					return false;
				}
			}
			wl.shrink_left(j);
		}
		// third: general constraint BCP
		if (wl.right_size() != 0) {
			WatchList::right_iterator it, end, j = wl.right_begin(); 
			Constraint::PropResult res;
			for (it = wl.right_begin(), end = wl.right_end(); it != end; ) {
				GenericWatch& w = *it++;
				res = w.propagate(*this, p);
				if (res.keepWatch) {
					*j++ = w;
				}
				if (!res.ok) {
					wl.shrink_right(std::copy(it, end, j));
					return false;
				}
			}
			wl.shrink_right(j);
		}
	}
	return DL || assign_.front == units_ || propagateUnits();
}

bool Solver::test(Literal p, PostPropagator* c) {
	assert(value(p.var()) == value_free && !hasConflict());
	assume(p); --stats.choices;
	freezeLevel(decisionLevel()); // can't split-off this level
	if (propagateUntil(c)) {
		assert(decision(decisionLevel()) == p);
		if (c) c->undoLevel(*this);
		undoUntil(decisionLevel()-1);
		return true;
	}
	unfreezeLevel(decisionLevel());
	assert(decision(decisionLevel()) == p);
	assign_.qReset();
	post_.reset();
	return false;
}

bool Solver::resolveConflict() {
	assert(hasConflict());
	if (decisionLevel() > rootLevel_) {
		if (decisionLevel() != btLevel_ && strategy_.search != SolverStrategies::no_learning) {
			uint32 uipLevel = analyzeConflict();
			stats.updateJumps(decisionLevel(), uipLevel, btLevel_, ccInfo_.lbd());
			undoUntil( uipLevel );
			return ClauseCreator::create(*this, cc_, ccInfo_);
		}
		else {
			return backtrack();
		}
	}
	return false;
}

bool Solver::backtrack() {
	Literal lastChoiceInverted;
	do {
		if (decisionLevel() == rootLevel_) return false;
		lastChoiceInverted = ~decision(decisionLevel());
		btLevel_ = decisionLevel() - 1;
		undoUntil(btLevel_, true);
	} while (hasConflict() || !force(lastChoiceInverted, 0));
	return true;
}

bool ImpliedList::assign(Solver& s) {
	assert(front <= lits.size());
	bool ok             = !s.hasConflict();
	const uint32 DL     = s.decisionLevel();
	VecType::iterator j = lits.begin() + front;
	for (VecType::iterator it = j, end = lits.end(); it != end; ++it) {
		if (it->level <= DL && ok) {
			ok = s.force(it->lit, it->ante.ante(), it->ante.data());
		}
		if (it->level < DL) { *j++ = *it; }
	}
	lits.erase(j, lits.end());
	level = DL * uint32(!lits.empty());
	front = level > s.rootLevel() ? front  : lits.size();
	return ok;
}

void Solver::undoUntil(uint32 level) {
	assert(btLevel_ >= rootLevel_);
	level      = std::max( level, btLevel_ );
	if (level >= decisionLevel()) return;
	uint32 num = decisionLevel() - level;
	bool sp    = strategy_.saveProgress > 0 && ((uint32)strategy_.saveProgress) <= num;
	bool ok    = conflict_.empty() && levels_.back().freeze == 0;
	conflict_.clear();
	heuristic_->undoUntil( *this, levels_[level].trailPos);
	undoLevel(sp && ok);
	while (--num) { undoLevel(sp); }
	assert(level == decisionLevel());
}

uint32 Solver::undoUntil(uint32 level, bool popBt) {
	if (popBt && backtrackLevel() > level && !shared_->project(decision(backtrackLevel()).var())) {
		setBacktrackLevel(level);
	}
	undoUntil(level);
	if (impliedLits_.active(level = decisionLevel())) {
		impliedLits_.assign(*this);
	}
	return level;
}

uint32 Solver::estimateBCP(const Literal& p, int rd) const {
	if (value(p.var()) != value_free) return 0;
	LitVec::size_type first = assign_.assigned();
	LitVec::size_type i     = first;
	Solver& self            = const_cast<Solver&>(*this);
	self.assign_.setValue(p.var(), trueValue(p));
	self.assign_.trail.push_back(p);
	const ShortImplicationsGraph& btig = shared_->shortImplications();
	do {
		Literal x = assign_.trail[i++];  
		if (!btig.propagateBin(self.assign_, x, 0)) {
			break;
		}
	} while (i < assign_.assigned() && rd-- != 0);
	i = assign_.assigned()-first;
	while (self.assign_.assigned() != first) {
		self.assign_.undoLast();
	}
	return (uint32)i;
}

uint32 Solver::inDegree(WeightLitVec& out) {
	if (decisionLevel() == 0) { return 1; }
	out.reserve((numAssignedVars()-levelStart(1))/10);
	uint32 maxIn  = 1;
	uint32 i      = assign_.trail.size(), stop = levelStart(1);
	for (Antecedent xAnte; i-- != stop; ) {
		Literal x     = assign_.trail[i];
		uint32  xLev  = assign_.level(x.var());
		xAnte         = assign_.reason(x.var());
		uint32  xIn   = 0;
		if (!xAnte.isNull() && xAnte.type() != Antecedent::binary_constraint) {
			conflict_.clear();
			xAnte.reason(*this, x, conflict_);
			for (LitVec::const_iterator it = conflict_.begin(); it != conflict_.end(); ++it) {
				xIn += level(it->var()) != xLev;
			}
			if (xIn) {
				out.push_back(WeightLiteral(x, xIn));
				maxIn     = std::max(xIn, maxIn);
			}
		}
	}
	return maxIn;
}
/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Private helper functions
////////////////////////////////////////////////////////////////////////////////////////
Solver::ConstraintDB* Solver::allocUndo(Constraint* c) {
	if (undoHead_ == 0) {
		return new ConstraintDB(1, c);
	}
	assert(undoHead_->size() == 1);
	ConstraintDB* r = undoHead_;
	undoHead_ = (ConstraintDB*)undoHead_->front();
	r->clear();
	r->push_back(c);
	return r;
}
void Solver::undoFree(ConstraintDB* x) {
	// maintain a single-linked list of undo lists
	x->clear();
	x->push_back((Constraint*)undoHead_);
	undoHead_ = x;
}
// removes the current decision level
void Solver::undoLevel(bool sp) {
	assert(decisionLevel() != 0 && levels_.back().trailPos != assign_.trail.size() && "Decision Level must not be empty");
	assign_.undoTrail(levels_.back().trailPos, sp);
	if (levels_.back().undo) {
		const ConstraintDB& undoList = *levels_.back().undo;
		for (ConstraintDB::size_type i = 0, end = undoList.size(); i != end; ++i) {
			undoList[i]->undoLevel(*this);
		}
		undoFree(levels_.back().undo);
	}
	levels_.pop_back();
}

inline ClauseHead* clause(const Antecedent& ante) {
	return ante.isNull() || ante.type() != Antecedent::generic_constraint ? 0 : ante.constraint()->clause();
}

// computes the First-UIP clause and stores it in cc_, where cc_[0] is the asserting literal (inverted UIP)
// and cc_[1] is a literal from the asserting level (if > 0)
// RETURN: asserting level of the derived conflict clause
uint32 Solver::analyzeConflict() {
	// must be called here, because we unassign vars during analyzeConflict
	heuristic_->undoUntil( *this, levels_.back().trailPos );
	uint32 onLevel  = 0;        // number of literals from the current DL in resolvent
	uint32 resSize  = 0;        // size of current resolvent
	Literal p;                  // literal to be resolved out next
	cc_.assign(1, p);           // will later be replaced with asserting literal
	Antecedent lhs, rhs, last;  // resolve operands
	const bool doOtfs = strategy_.otfs > 0;
	for (bumpAct_.clear();;) {
		uint32 lhsSize = resSize;
		uint32 rhsSize = 0;
		heuristic_->updateReason(*this, conflict_, p);
		for (LitVec::size_type i = 0; i != conflict_.size(); ++i) {
			Literal& q = conflict_[i];
			uint32 cl  = level(q.var());
			rhsSize   += (cl != 0);
			if (!seen(q.var())) {
				++resSize;
				assert(isTrue(q) && "Invalid literal in reason set!");
				assert(cl > 0 && "Top-Level implication not marked!");
				markSeen(q.var());
				if (cl == decisionLevel()) {
					++onLevel;
				}
				else {
					cc_.push_back(~q);
					markLevel(cl);
				}
			}
		}
		if (resSize != lhsSize) { lhs = 0; }
		if (rhsSize != resSize) { rhs = 0; }
		if (doOtfs && (!rhs.isNull() || !lhs.isNull())) {
			// resolvent subsumes rhs and possibly also lhs
			otfs(lhs, rhs, p, onLevel == 1);
		}
		assert(onLevel > 0 && "CONFLICT MUST BE ANALYZED ON CONFLICT LEVEL!");
		// search for the last assigned literal that needs to be analyzed...
		while (!seen(assign_.last().var())) {
			assign_.undoLast();
		}
		p   = assign_.last();
		rhs = reason(p);
		clearSeen(p.var());
		if (--onLevel == 0) { 
			break; 
		}
		--resSize; // p will be resolved out next
		last = rhs;
		reason(p, conflict_);
	}
	cc_[0] = ~p; // store the 1-UIP
	assert(decisionLevel() == level(cc_[0].var()));
	ClauseHead* lastRes = 0;
	if (strategy_.otfs > 1 || !lhs.isNull()) {
		if (!lhs.isNull()) { 
			lastRes = clause(lhs);
		}
		else if (cc_.size() <= (conflict_.size()+1)) {
			lastRes = clause(last);
		}
	}
	if (strategy_.bumpVarAct && reason(p).learnt()) {
		bumpAct_.push_back(WeightLiteral(p, static_cast<LearntConstraint*>(reason(p).constraint())->activity().lbd()));
	}
	ccInfo_ = ClauseInfo(Constraint_t::learnt_conflict);
	return simplifyConflictClause(cc_, ccInfo_, lastRes);
}

void Solver::otfs(Antecedent& lhs, const Antecedent& rhs, Literal p, bool final) {
	ClauseHead* cLhs = clause(lhs), *cRhs = clause(rhs);
	ClauseHead::BoolPair x;
	if (cLhs) {
		x = cLhs->strengthen(*this, ~p, !final);
		if (!x.first || x.second) {
			cLhs = !x.first ? 0 : otfsRemove(cLhs, 0);
		}
	}
	lhs = cLhs;
	if (cRhs) {
		x = cRhs->strengthen(*this, p, !final);
		if (!x.first || (x.second && otfsRemove(cRhs, 0) == 0)) {
			if (x.first && reason(p) == cRhs) { setReason(p, 0); }
			cRhs = 0;
		}
		if (cLhs && cRhs) {
			// lhs and rhs are now equal - only one of them is needed
			if (!cLhs->learnt()) {
				std::swap(cLhs, cRhs);
			}
			otfsRemove(cLhs, 0);
		}
		lhs = cRhs;
	}
}

ClauseHead* Solver::otfsRemove(ClauseHead* c, const LitVec* newC) {
	bool remStatic = !newC || (newC->size() <= 3 && shared_->allowImplicit(Constraint_t::learnt_conflict));
	if (c->learnt() || remStatic) {
		ConstraintDB& db = (c->learnt() ? learnts_ : constraints_);
		ConstraintDB::iterator it;
		if ((it = std::find(db.begin(), db.end(), c)) != db.end()) {
			db.erase(it);
			c->destroy(this, true);
			c = 0;
		}
	}
	return c;
}

// minimizes the conflict clause in cc w.r.t selected strategies.
// PRE:
//  - cc is a valid conflict clause and cc[0] is the UIP-literal
//  - all literals in cc except cc[0] are marked
//  - all decision levels of literals in cc are marked
//  - rhs is 0 or a clause that might be subsumed by cc
// RETURN: finalizeConflictClause(cc, info)
uint32 Solver::simplifyConflictClause(LitVec& cc, ClauseInfo& info, ClauseHead* rhs) {
	// 1. remove redundant literals from conflict clause
	temp_.clear();
	if (!ccMin_ && strategy_.strRecursive) { ccMin_ = new CCMinRecursive; }
	uint32 onAssert = ccMinimize(cc, temp_, strategy_.ccMinAntes, ccMin_);
	uint32 jl       = cc.size() > 1 ? level(cc[1].var()) : 0;
	// clear seen flags of removed literals - keep levels marked
	for (LitVec::size_type x = 0, stop = temp_.size(); x != stop; ++x) {
		clearSeen(temp_[x].var());
	}
	// 2. check for inverse arcs
	if (onAssert == 1 && strategy_.reverseArcs > 0) {
		uint32 maxN = (uint32)strategy_.reverseArcs;
		if      (maxN > 2) maxN = UINT32_MAX;
		else if (maxN > 1) maxN = static_cast<uint32>(cc.size() / 2);
		markSeen(cc[0].var());
		Antecedent ante = ccHasReverseArc(cc[1], jl, maxN);
		if (!ante.isNull()) {
			// resolve with inverse arc
			conflict_.clear();
			ante.reason(*this, ~cc[1], conflict_);
			ccResolve(cc, 1, conflict_);
		}
		clearSeen(cc[0].var());
	}
	// 3. check if final clause subsumes rhs
	if (rhs) {
		conflict_.clear();
		rhs->toLits(conflict_);
		uint32 open   = (uint32)cc.size();
		markSeen(cc[0].var());
		for (LitVec::const_iterator it = conflict_.begin(), end = conflict_.end(); it != end && open; ++it) {
			// NOTE: at this point the DB might not be fully simplified,
			//       e.g. because of mt or lookahead, hence we must explicitly
			//       check for literals assigned on DL 0
			open -= level(it->var()) > 0 && seen(it->var());
		}
		rhs = open ? 0 : otfsRemove(rhs, &cc); 
		if (rhs) { // rhs is subsumed by cc but could not be removed.
			// TODO: we could reuse rhs instead of learning cc
			//       but this would complicate the calling code.
			ClauseHead::BoolPair r(true, false);
			if (cc_.size() < conflict_.size()) {
				//     For now, we only try to strengthen rhs.
				for (LitVec::const_iterator it = conflict_.begin(), end = conflict_.end(); it != end && r.first; ++it) {
					if (!seen(it->var()) || level(it->var()) == 0) {
						r = rhs->strengthen(*this, *it, false);
					}
				}
				if (!r.first) { rhs = 0; }
			}
		}
		clearSeen(cc[0].var());
	}
	// 4. finalize
	jl = finalizeConflictClause(cc, info);
	// 5. bump vars implied by learnt constraints with small lbd 
	if (!bumpAct_.empty()) {
		WeightLitVec::iterator j = bumpAct_.begin();
		weight_t newLbd = info.lbd();
		for (WeightLitVec::iterator it = bumpAct_.begin(), end = bumpAct_.end(); it != end; ++it) {
			if (it->second < newLbd) {
				it->second = 1 + (it->second <= 2); 
				*j++ = *it;	
			}
		}
		bumpAct_.erase(j, bumpAct_.end());
		heuristic_->bump(*this, bumpAct_, 1.0);
	}
	bumpAct_.clear();
	// 6. clear level flags of redundant literals
	for (LitVec::size_type x = 0, stop = temp_.size(); x != stop; ++x) {
		unmarkLevel(level(temp_[x].var()));
	}
	return jl;
}

// conflict clause minimization
// PRE: 
//  - cc is an asserting clause and cc[0] is the asserting literal
//  - all literals in cc are marked as seen
//  -  if ccMin != 0, all decision levels of literals in cc are marked
// POST:
//  - redundant literals were added to removed
//  - if (cc.size() > 1): cc[1] is a literal from the asserting level
// RETURN
//  - the number of literals from the asserting level
uint32 Solver::ccMinimize(LitVec& cc, LitVec& removed, uint32 antes, CCMinRecursive* ccMin) {
	if (ccMin) { ccMin->init(numVars()+1); }
	// skip the asserting literal
	LitVec::size_type j = 1;
	uint32 assertLevel  = 0;
	uint32 assertPos    = 1;
	uint32 onAssert     = 0;
	uint32 varLevel     = 0;
	for (LitVec::size_type i = 1; i != cc.size(); ++i) { 
		if (antes == 0 || !ccRemovable(~cc[i], antes-1, ccMin)) {
			if ( (varLevel = level(cc[i].var())) > assertLevel ) {
				assertLevel = varLevel;
				assertPos   = static_cast<uint32>(j);
				onAssert    = 0;
			}
			onAssert += (varLevel == assertLevel);
			cc[j++] = cc[i];
		}
		else { 
			removed.push_back(cc[i]);
		}
	}
	cc.erase(cc.begin()+j, cc.end());
	if (assertPos != 1) {
		std::swap(cc[1], cc[assertPos]);
	}
	if (ccMin) { ccMin->clear(); }
	return onAssert;
}

// returns true if p is redundant in current conflict clause
bool Solver::ccRemovable(Literal p, uint32 antes, CCMinRecursive* ccMin) {
	const Antecedent& ante = reason(p);
	if (ante.isNull() || !(antes <= (uint32)ante.type())) {
		return false;
	}
	if (!ccMin) { return ante.minimize(*this, p, 0); }
	// recursive minimization
	LitVec& dfsStack = ccMin->dfsStack;
	assert(dfsStack.empty());
	CCMinRecursive::State dfsState = CCMinRecursive::state_removable;
	p.clearWatch();
	dfsStack.push_back(p);
	for (Literal x;; ) {
		x = dfsStack.back();
		dfsStack.pop_back();
		assert(!seen(x.var()) || x == p);
		if (x.watched()) {
			if (x == p) return dfsState == CCMinRecursive::state_removable;
			ccMin->markVisited(x, dfsState);
		}
		else if (dfsState != CCMinRecursive::state_poison) {
			CCMinRecursive::State temp = ccMin->state(x);
			if (temp == CCMinRecursive::state_open) {
				assert(value(x.var()) != value_free && hasLevel(level(x.var())));
				x.watch();
				dfsStack.push_back(x);
				const Antecedent& ante = reason(x);
				if (ante.isNull() || !(antes <= (uint32)ante.type()) || !ante.minimize(*this, x, ccMin)) {
					dfsState = CCMinRecursive::state_poison;
				}
			}
			else if (temp == CCMinRecursive::state_poison) {
				dfsState = temp;
			}
		}
	}
}

// checks whether there is a valid "inverse arc" for the given literal p that can be used
// to resolve p out of the current conflict clause
// PRE: 
//  - all literals in the current conflict clause are marked
//  - p is a literal of the current conflict clause and level(p) == maxLevel
// RETURN
//  - An antecedent that is an "inverse arc" for p or null if no such antecedent exists.
Antecedent Solver::ccHasReverseArc(const Literal p, uint32 maxLevel, uint32 maxNew) {
	assert(seen(p.var()) && isFalse(p) && level(p.var()) == maxLevel);
	const ShortImplicationsGraph& btig = shared_->shortImplications();
	Antecedent ante;
	if (btig.reverseArc(*this, p, maxLevel, ante)) { return ante; }
	WatchList& wl   = watches_[p.index()];
	WatchList::left_iterator it, end; 
	for (it = wl.left_begin(), end = wl.left_end();  it != end;  ++it) {
		if (it->head->isReverseReason(*this, ~p, maxLevel, maxNew)) {
			return it->head;
		}
	}
	return ante;
}

// removes cc[pos] by resolving cc with reason
void Solver::ccResolve(LitVec& cc, uint32 pos, const LitVec& reason) {
	heuristic_->updateReason(*this, reason, cc[pos]);
	Literal x;
	for (LitVec::size_type i = 0; i != reason.size(); ++i) {
		x = reason[i];
		assert(isTrue(x));
		if (!seen(x.var())) {
			markLevel(level(x.var()));
			cc.push_back(~x);
		}
	}
	clearSeen(cc[pos].var());
	unmarkLevel(level(cc[pos].var()));
	cc[pos] = cc.back();
	cc.pop_back();
}

// computes asserting level and lbd of cc and clears flags.
// POST:
//  - literals and decision levels in cc are no longer marked
//  - if cc.size() > 1: cc[1] is a literal from the asserting level
// RETURN: asserting level of conflict clause.
uint32 Solver::finalizeConflictClause(LitVec& cc, ClauseInfo& info) {
	// 2. clear flags and compute lbd
	uint32  lbd         = 1;
	uint32  onRoot      = 0;
	uint32  varLevel    = 0;
	uint32  assertLevel = 0;
	uint32  assertPos   = 1;
	Literal tagLit      = ~sharedContext()->tagLiteral();
	for (LitVec::size_type i = 1; i != cc.size(); ++i) {
		clearSeen(cc[i].var());
		if (cc[i] == tagLit) { info.setTagged(true); }
		if ( (varLevel = level(cc[i].var())) > assertLevel ) {
			assertLevel = varLevel;
			assertPos   = static_cast<uint32>(i);
		}
		if (hasLevel(varLevel)) {
			unmarkLevel(varLevel);
			lbd += (varLevel > rootLevel()) || (++onRoot == 1);
		}
	}
	if (assertPos != 1) { std::swap(cc[1], cc[assertPos]); }
	info.setLbd(lbd);
	info.setActivity(static_cast<uint32>(1+stats.restarts));
	return assertLevel;
}

// (inefficient) default implementation 
bool Constraint::minimize(Solver& s, Literal p, CCMinRecursive* rec) {
	LitVec temp;
	reason(s, p, temp);
	for (LitVec::size_type i = 0; i != temp.size(); ++i) {
		if (!s.ccMinimize(temp[i], rec)) {
			return false;
		}
	}
	return true;
}

// Selects next branching literal
bool Solver::decideNextBranch(double f) { 
	if (f <= 0.0 || strategy_.rng.drand() >= f || numFreeVars() == 0) {
		return heuristic_->select(*this);
	}
	// select randomly
	Literal choice;
	uint32 maxVar = numVars() + 1; 
	for (uint32 v = strategy_.rng.irand(maxVar);;) {
		if (value(v) == value_free) {
			choice    = DecisionHeuristic::selectLiteral(*this, v, 0);
			break;
		}
		if (++v == maxVar) { v = 1; }
	}
	return assume(choice);
}
// Removes up to remFrac% of the learnt nogoods but
// keeps those that are locked or are highly active.
Solver::DBInfo Solver::reduceLearnts(float remFrac, const ReduceStrategy& rs) {
	uint32 oldS = numLearntConstraints();
	uint32 remM = static_cast<uint32>(oldS * std::max(0.0f, remFrac));
	DBInfo r    = {0,0,0};
	CmpScore cmp(learnts_, (ReduceStrategy::Score)rs.score, rs.glue);
	if (remM >= oldS || !remM || rs.algo == ReduceStrategy::reduce_sort) {
		r = reduceSortInPlace(remM, cmp, false);
	}
	else if (rs.algo == ReduceStrategy::reduce_stable) { r = reduceSort(remM, cmp);  }
	else if (rs.algo == ReduceStrategy::reduce_heap)   { r = reduceSortInPlace(remM, cmp, true);}
	else                                               { r = reduceLinear(remM, cmp); }
	stats.addDeleted(oldS - r.size);
	shrinkVecTo(learnts_, r.size);
	return r;
}

// Removes up to maxR of the learnt nogoods.
// Keeps those that are locked or have a high activity and
// does not reorder learnts_.
Solver::DBInfo Solver::reduceLinear(uint32 maxR, const CmpScore& sc) {
	// compute average activity
	uint64 scoreSum  = 0;
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		scoreSum += sc.score(static_cast<LearntConstraint*>(learnts_[i])->activity());
	}
	double avgAct    = (scoreSum / (double) numLearntConstraints());
	// constraints with socre > 1.5 times the average are "active"
	double scoreThresh = avgAct * 1.5;
	double scoreMax    = (double)sc.score(Activity(Activity::MAX_ACT, 1));
	if (scoreThresh > scoreMax) {
		scoreThresh = (scoreMax + (scoreSum / (double) numLearntConstraints())) / 2.0;
	}
	// remove up to maxR constraints but keep "active" and locked once
	const uint32 glue = sc.glue;
	DBInfo       res  = {0,0,0};
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		LearntConstraint* c = static_cast<LearntConstraint*>(learnts_[i]);
		Activity a          = c->activity();
		bool isLocked       = c->locked(*this);
		bool isGlue         = (sc.score(a) > scoreThresh || a.lbd() <= glue);
		if (maxR == 0 || isLocked || isGlue) {
			res.pinned             += isGlue;
			res.locked             += isLocked;
			learnts_[res.size++]    = c;
			c->decreaseActivity();
		}
		else {
			--maxR;
			c->destroy(this, true);
		}
	}
	return res;
}

// Sorts learnt constraints by score and removes the
// maxR constraints with the lowest score without
// reordering learnts_.
Solver::DBInfo Solver::reduceSort(uint32 maxR, const CmpScore& sc) {
	typedef PodVector<CmpScore::ViewPair>::type HeapType;
	maxR              = std::min(maxR, (uint32)learnts_.size());
	const uint32 glue = sc.glue;
	DBInfo       res  = {0,0,0};
	HeapType     heap;
	heap.reserve(maxR);
	bool isGlue, isLocked;
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		LearntConstraint* c = static_cast<LearntConstraint*>(learnts_[i]);
		CmpScore::ViewPair vp(i, c->activity());
		res.pinned += (isGlue   = vp.second.lbd() <= glue);
		res.locked += (isLocked = c->locked(*this));
		if (!isLocked && !isGlue) { 
			if (maxR) { // populate heap with first maxR constraints
				heap.push_back(vp); 
				if (--maxR == 0) { std::make_heap(heap.begin(), heap.end(), sc); } 
			}
			else if (sc(vp, heap[0])) { // replace max element in heap
				std::pop_heap(heap.begin(), heap.end(), sc);
				heap.back() = vp;
				std::push_heap(heap.begin(), heap.end(), sc);
			}
		}
	}
	// Remove all constraints in heap - those are "inactive".
	for (HeapType::const_iterator it = heap.begin(), end = heap.end(); it != end; ++it) {
		learnts_[it->first]->destroy(this, true);
		learnts_[it->first] = 0;
	}
	// Cleanup db and decrease activity of remaining constraints.
	uint32 j = 0;
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		if (LearntConstraint* c = static_cast<LearntConstraint*>(learnts_[i])) {
			c->decreaseActivity();
			learnts_[j++] = c;
		}
	}
	res.size = j;
	return res;
}

// Sorts the learnt db by score and removes the first 
// maxR constraints (those with the lowest score). 
Solver::DBInfo Solver::reduceSortInPlace(uint32 maxR, const CmpScore& sc, bool partial) {
	maxR                        = std::min(maxR, (uint32)learnts_.size());
	ConstraintDB::iterator nEnd = learnts_.begin();
	const uint32 glue           = sc.glue;
	DBInfo res                  = {0,0,0};
	bool isGlue, isLocked;
	if (!partial) {
		// sort whole db and remove first maxR constraints
		if (maxR && maxR != learnts_.size()) std::stable_sort(learnts_.begin(), learnts_.end(), sc);
		for (ConstraintDB::iterator it = learnts_.begin(), end = learnts_.end(); it != end; ++it) {
			LearntConstraint* c = static_cast<LearntConstraint*>(*it);
			res.pinned         += (isGlue = c->activity().lbd() <= glue); 
			res.locked         += (isLocked = c->locked(*this));
			if (!maxR || isLocked || isGlue) { c->decreaseActivity(); *nEnd++ = c; }
			else                             { c->destroy(this, true); --maxR; }
		}
	}
	else {
		// reorder db s.th first maxR constraints are a heap of most inactive
		ConstraintDB::iterator hBeg = learnts_.begin();
		ConstraintDB::iterator hEnd = learnts_.begin();
		for (ConstraintDB::iterator it = learnts_.begin(), end = learnts_.end(); it != end; ++it) {
			LearntConstraint* c = static_cast<LearntConstraint*>(*it);
			res.pinned         += (isGlue = c->activity().lbd() <= glue); 
			res.locked         += (isLocked = c->locked(*this));
			if      (isLocked || isGlue) { continue; }
			else if (maxR)               {
				*it     = *hEnd;
				*hEnd++ = c;
				if (--maxR == 0) { std::make_heap(hBeg, hEnd, sc); }
			}
			else if (sc(c, learnts_[0])) {
				std::pop_heap(hBeg, hEnd, sc);
				*it      = *(hEnd-1);
				*(hEnd-1)= c;
				std::push_heap(hBeg, hEnd, sc);
			}
		}
		// remove all constraints in heap
		for (ConstraintDB::iterator it = hBeg; it != hEnd; ++it) {
			static_cast<LearntConstraint*>(*it)->destroy(this, true);
		}
		// copy remaining constraints down
		for (ConstraintDB::iterator it = hEnd, end = learnts_.end(); it != end; ++it) {
			LearntConstraint* c = static_cast<LearntConstraint*>(*it);
			c->decreaseActivity();
			*nEnd++ = c;
		}
	}
	res.size = static_cast<uint32>(std::distance(learnts_.begin(), nEnd));
	return res;	
}

uint32 Solver::countLevels(const Literal* first, const Literal* last, uint32 maxLevel) {
	if (maxLevel < 2)    { return uint32(maxLevel && first != last); }
	if (++lbdTime_ != 0) { lbdStamp_.resize(levels_.size()+1, lbdTime_-1); }
	else                 { lbdStamp_.assign(levels_.size()+1, lbdTime_); lbdTime_ = 1; }
	lbdStamp_[0] = lbdTime_;
	uint32 levels= 0;
	for (uint32 lev; first != last; ++first) {
		lev = level(first->var());
		if (lbdStamp_[lev] != lbdTime_) {
			lbdStamp_[lev] = lbdTime_;
			if (++levels == maxLevel) { break; }
		}
	}
	return levels;
}

uint64 Solver::updateBranch(uint32 cfl) {
	int32 dl = (int32)decisionLevel(), xl = static_cast<int32>(cflStamp_.size())-1;
	if      (xl > dl) { do { cfl += cflStamp_.back(); cflStamp_.pop_back(); } while (--xl != dl); }
	else if (dl > xl) { cflStamp_.insert(cflStamp_.end(), dl - xl, 0); }
	return cflStamp_.back() += cfl;
}
/////////////////////////////////////////////////////////////////////////////////////////
// The basic DPLL-like search-function
/////////////////////////////////////////////////////////////////////////////////////////
ValueRep Solver::search(SearchLimits& limit, double rf) {
	uint64 local = limit.local != UINT64_MAX ? limit.local : 0;
	rf           = std::max(0.0, std::min(1.0, rf));
	if (local && decisionLevel() == rootLevel()) { cflStamp_.assign(decisionLevel()+1, 0); }
	do {
		for (uint32 conflicts = hasConflict() || !propagate() || !simplify();;) {
			if (conflicts) {
				for (conflicts = 1; resolveConflict() && !propagate(); ) { ++conflicts; }
				limit.conflicts -= conflicts < limit.conflicts ? conflicts : limit.conflicts;
				if (local && updateBranch(conflicts) >= local)                  { limit.local = 0; }
				if (hasConflict() || (decisionLevel() == 0 && !simplify()))     { return value_false; }
				if ((limit.reached(stats) || learntLimit(limit))&&numFreeVars()){ return value_free;  }
			}
			if (decideNextBranch(rf)) { conflicts = !propagate(); }
			else                      { break; }
		}
	} while (!post_.isModel(*this));
	if (satPrepro()) { temp_.clear(); satPrepro()->extendModel(assign_, temp_); }
	return value_true;
}

bool Solver::nextSymModel() {
	assert(numFreeVars() == 0);
	if (satPrepro() && !temp_.empty()) {
		satPrepro()->extendModel(assign_, temp_);
		return true;
	}
	return false;
}
}
