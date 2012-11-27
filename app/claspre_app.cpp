// 
// Copyright (c) 2012, Benjamin Kaufmann
// 
// This file is part of Claspre. See http://potassco.sourceforge.net/
// 
// Claspre is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Claspre is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Claspre; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
#include "claspre_app.h"
#include "alarm.h"
#include <program_opts/composite_value_parser.h>  // pair and vector
#include <program_opts/typed_value.h>
#include <clasp/satelite.h>
#include <clasp/minimize_constraint.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>
#include <math.h>
#if !defined(_WIN32)
#include <unistd.h> // for _exit
#endif
using namespace Clasp;
namespace Claspre {
/////////////////////////////////////////////////////////////////////////////////////////
// Output
/////////////////////////////////////////////////////////////////////////////////////////

// clasp is done - write result
void Output::onExit(const Solver& s, const Result& r) {

	printDynamic(s);

	if (s.sharedContext()->enumerator()->minimize()) {
		printf(",\n");
		printf(" \"Optimization\": [\n");
		printf("  [\"Avg_Improvement\", %.4f],\n", ostats.avg_impr);
		double var_impr = s.stats.models > 1 ? ostats.var_impr / (s.stats.models - 1): 0;
		printf("  [\"Stdev_Improvement\", %.4f],\n", sqrt(var_impr));
		printf("  [\"Var_Coeff_Improvement\", %.4f],\n", sqrt(var_impr) / ostats.avg_impr);
		printf("  [\"Avg_Ratio_Improvement\", %.4f],\n", ostats.avg_ratio_impr);
		double var_ration_impr = s.stats.models > 1 ? ostats.var_ratio_impr / (s.stats.models - 1): 0;
		printf("  [\"Stdev_Ratio_Improvement\", %.4f],\n", sqrt(var_ration_impr));
		printf("  [\"Var_Coeff_Ratio_Improvement\", %.4f],\n", sqrt(var_ration_impr) / ostats.avg_ratio_impr);
		printf("  [\"Ratio_Worst_Best\", %.4f]\n", static_cast<double>(ostats.first_quality) / ostats.last_quality );
		printf(" ]\n");
	}
	printf("}\n");
//	printf("Times:\n");
//	printf("  Total: %.3f\n", r.totalTime);
//	printf("  Solve: %.3f\n", r.solveTime);
}

// state change in clasp
void Output::onState(bool enter) {
	// <TBD>
	// add your code here
	//TODO: Nothing to do here?
	// <TBD>
	if (facade_.state() == ClaspFacade::state_read && enter) {
		//printf("Reading from %s\n", opts_.input[0].c_str());
	}
	if (facade_.state() == ClaspFacade::state_preprocess && enter) {
		//printf("Preprocessing\n");
	}

}

// preprocessing done
void Output::onProgramPrepared(const Solver& s) {
	//printf("Preprocessing done!\n");
	printf("{\n");
	const PreproStats* logicProgram = facade_.api() ? &facade_.api()->stats : 0;
	if (logicProgram) {
		// ASP stats
		uint32 tight = logicProgram->sccs == 0 || logicProgram->sccs == PrgNode::noScc;
		uint32 sccs  = logicProgram->sccs != PrgNode::noScc ? logicProgram->sccs : 0;

		printf(" \"After_Preprocessing\": [\n");
		printf("  [\"Tight\", %u],\n",  tight);
		printf("  [\"Problem_Variables\", %u],\n",  s.numVars());
		printf("  [\"Free_Problem_Variables\", %u],\n",  s.numFreeVars());
		printf("  [\"Assigned_Problem_Variables\", %u],\n",  s.numAssignedVars());
		printf("  [\"Constraints\", %u],\n",  s.numConstraints());
		printf("  [\"Constraints/Vars\", %.4f],\n",  static_cast<double>(s.numConstraints()) / s.numVars());
		printf("  [\"Created_Bodies\", %u],\n", logicProgram->bodies);
		printf("  [\"Program_Atoms\", %u],\n", logicProgram->atoms);
		printf("  [\"SCCS\", %u],\n", logicProgram->sccs);
		printf("  [\"Non-trivial_SCCS\", %u],\n", sccs);
		printf("  [\"Nodes_in_Positive_BADG\", %u],\n", logicProgram->ufsNodes);
		printf("  [\"Rules\", %u],\n", logicProgram->rules[0]);
		printf("  [\"Normal_Rules\", %u],\n", logicProgram->rules[1]);
		printf("  [\"Cardinality_Rules\", %u],\n", logicProgram->rules[2]);
		printf("  [\"Choice_Rules\", %u],\n", logicProgram->rules[3]);
		printf("  [\"Weight_Rules\", %u],\n", logicProgram->rules[5]);
		printf("  [\"Optimization_Rules\", %u],\n", logicProgram->rules[6]);
		printf("  [\"Frac_Normal_Rules\", %.4f],\n", static_cast<double>(logicProgram->rules[1]) / logicProgram->rules[0]);
		printf("  [\"Frac_Cardinality_Rules\", %.4f],\n", static_cast<double>(logicProgram->rules[2]) / logicProgram->rules[0]);
		printf("  [\"Frac_Choice_Rules\", %.4f],\n", static_cast<double>(logicProgram->rules[3]) / logicProgram->rules[0]);
		printf("  [\"Frac_Weight_Rules\", %.4f],\n", static_cast<double>(logicProgram->rules[5]) / logicProgram->rules[0]);
		printf("  [\"Frac_Optimization_Rules\", %.4f],\n", static_cast<double>(logicProgram->rules[6]) / logicProgram->rules[0]);
		printf("  [\"Equivalences\", %u],\n", logicProgram->sumEqs());
		printf("  [\"Atom-Atom_Equivalences\", %u],\n", logicProgram->eqs[0]);
		printf("  [\"Body-Body_Equivalences\", %u],\n", logicProgram->eqs[1]);
		printf("  [\"Other_Equivalences\", %u],\n", logicProgram->eqs[2]);
		if (logicProgram->sumEqs() > 0) {
			printf("  [\"Frac_Atom-Atom_Equivalences\", %.4f],\n",  static_cast<double>(logicProgram->eqs[0]) /logicProgram->sumEqs());
			printf("  [\"Frac_Body-Body_Equivalences\", %.4f],\n",  static_cast<double>(logicProgram->eqs[1]) /logicProgram->sumEqs());
			printf("  [\"Frac_Other_Equivalences\", %.4f],\n",  static_cast<double>(logicProgram->eqs[2]) /logicProgram->sumEqs());
		}
		else {
			printf("  [\"Frac_Atom-Atom_Equivalences\", 0],\n");
			printf("  [\"Frac_Body-Body_Equivalences\", 0],\n");
			printf("  [\"Frac_Other_Equivalences\", 0],\n");
		}


		// Shared Context Stats
		printf("  [\"Binary_Constraints\", %u],\n",  s.sharedContext()->numBinary());
		printf("  [\"Ternary_Constraints\", %u],\n",  s.sharedContext()->numTernary());
		uint32 other_const = s.numConstraints() - s.sharedContext()->numTernary();
		printf("  [\"Other_Constraints\", %u],\n",  other_const);
		printf("  [\"Frac_Binary_Constraints\", %.4f],\n",  static_cast<double>(s.sharedContext()->numBinary()) / s.sharedContext()->numConstraints());
		printf("  [\"Frac_Ternary_Constraints\", %.4f],\n",  static_cast<double>(s.sharedContext()->numTernary()) / s.sharedContext()->numConstraints());
		printf("  [\"Frac_Other_Constraints\", %.4f]\n",  static_cast<double>(other_const) / s.sharedContext()->numConstraints());

		printf(" ]\n");

	}
}

void Output::onRestart(const Solver& s, uint64 conflictLimit, uint32 learntLimit) {
	printDynamic(s);
	printf("  ,\n");
}

void Output::printDynamic(const Solver& s) {
	const SolveStats& stats = s.stats;
	//printf("(Re)-Start %" PRIu64 " with limits (%" PRIu64 ", %u)\n", stats.restarts, conflictLimit, learntLimit);
	if (stats.restarts != 0) {
		//Core Stats
		printf(" \"Restart-%" PRIu64 "\" :[\n", stats.restarts);
		printf("  [\"Models\" , %" PRIu64 "],\n", stats.models);	//for optimization probs
		printf("  [\"Choices\" , %" PRIu64 "],\n", stats.choices);
		printf("  [\"Analyed_Conflicts\" , %" PRIu64 "],\n", stats.analyzed);
		printf("  [\"Conflicts/Choices\" , %.4f],\n", static_cast<double>(stats.analyzed) / stats.choices);
		printf("  [\"Avg_Conflict_Levels\" , %.4f],\n", stats.avgCfl());
		printf("  [\"Avg_LBD_Levels\" , %.4f],\n", stats.avgLbd());
		printf("  [\"Learnt_from_Conflict\" , %" PRIu64 "],\n", stats.learnts[0]);
		printf("  [\"Learnt_from_Loop\" , %" PRIu64 "],\n", stats.learnts[1]); //unfounded set checking
		printf("  [\"Learnt_from_Other\" , %" PRIu64 "],\n", stats.learnts[2]);
		uint64 sum_learnts = stats.learnts[0] + stats.learnts[1] + stats.learnts[2];
		printf("  [\"Frac_Learnt_from_Conflict\" , %.4f],\n", static_cast<double>(stats.learnts[0]) / sum_learnts);
		printf("  [\"Frac_Learnt_from_Loop\" , %.4f],\n", static_cast<double>(stats.learnts[1]) / sum_learnts); //unfounded set checking
		printf("  [\"Frac_Learnt_from_Other\" , %.4f],\n", static_cast<double>(stats.learnts[2]) / sum_learnts);
		printf("  [\"Literals_in_Conflict_Nogoods\" , %" PRIu64 "],\n", stats.lits[0]);
		printf("  [\"Literals_in_Loop_Nogoods\" , %" PRIu64 "],\n", stats.lits[1]); //unfounded set checking
		printf("  [\"Literals_in_Other_Nogoods\" , %" PRIu64 "],\n", stats.lits[2]);
		uint64 sum_lits = stats.lits[0] + stats.lits[1] +stats.lits[2];
		printf("  [\"Literals_in_Conflict_Nogoods\" , %.4f],\n", static_cast<double>(stats.lits[0]) / sum_lits);
		printf("  [\"Literals_in_Loop_Nogoods\" , %.4f],\n", static_cast<double>(stats.lits[1]) / sum_lits); //unfounded set checking
		printf("  [\"Literals_in_Other_Nogoods\" , %.4f],\n", static_cast<double>(stats.lits[2]) / sum_lits);
		printf("  [\"Removed_Nogood\" , %" PRIu64 "],\n", stats.deleted);
		printf("  [\"Learnt_Binary\" , %u],\n", stats.binary);
		printf("  [\"Learnt_Ternary\" , %u],\n", stats.ternary);
		printf("  [\"Frac_Removed_Nogood\" , %.4f],\n", static_cast<double>(stats.deleted) / sum_learnts);
		printf("  [\"Frac_Learnt_Binary\" , %.4f],\n", static_cast<double>(stats.binary) / sum_learnts);
		printf("  [\"Frac_Learnt_Ternary\" , %.4f],\n", static_cast<double>(stats.ternary) / sum_learnts);

		// Jump Stats
		const JumpStats* jstats = stats.jumps;
		printf("  [\"Decision_Literals_Models\" , %" PRIu64 "],\n", jstats->modLits); //for optimization problems
		printf("  [\"Skipped_Levels_while_Backjumping\" , %" PRIu64 "],\n", jstats->jumpSum);
		printf("  [\"Avg_Skipped_Levels_while_Backjumping\" , %.4f],\n", static_cast<double>(jstats->jumpSum) / stats.analyzed);
		printf("  [\"Longest_Backjumping\" , %u],\n", jstats->maxJump);
		// Averages Stats
		const SumQueue* sstats = stats.queue;
		printf("  [\"Running_Avg_Conflictlevel\" , %.4f],\n", sstats->avgCfl());
		printf("  [\"Running_Avg_LBD\" , %.4f]\n", sstats->avgLbd());
		printf(" ]\n");
	}
}

void Output::onDeletion(const Solver& s, uint64 conflictLimit, uint32 learntLimit) {
	printf("Deletion: next limits (%" PRIu64 ", %u)\n", conflictLimit, learntLimit);
	// <TBD>
}

// called on each model
void Output::onModel(const Solver& s, const Clasp::Enumerator& en) {
	/*printf("Model %" PRIu64 " on DL: %u\n", en.enumerated, s.decisionLevel());
	if (en.minimize()) {
		const SharedMinimizeData::SumVec& opt = en.minimize()->optimum()->opt;
		printf("Opt: %" PRId64, opt[0]);
		for (uint32 i = 1, end = en.minimize()->numRules(); i != end; ++i) {
			printf(" %" PRId64, opt[i]);
		}
		printf("\n");
	}*/
	if (en.minimize()) {
		const SharedMinimizeData::SumVec& opt = en.minimize()->optimum()->opt;
		uint64 index = en.enumerated;
		uint64 quality = opt[0];
		if (index == 1) {
			ostats.avg_impr  = 0.0;
			ostats.avg_ratio_impr = 0.0;
			ostats.var_impr = 0.0;
			ostats.first_quality = quality;
		}
		if (index > 1) {
			uint64 impr = ostats.last_quality - quality;
			double delta = impr - ostats.avg_impr;
			ostats.avg_impr = ostats.avg_impr + delta / static_cast<double>(index - 1.0);
			ostats.var_impr = ostats.var_impr + delta * (impr - ostats.avg_impr);

			double ratio_impr = static_cast<double>(ostats.last_quality) / quality;
			double delta_ratio = ratio_impr - ostats.avg_ratio_impr;
			ostats.avg_ratio_impr = ostats.avg_ratio_impr + delta_ratio / static_cast<double>(index - 1.0);
			ostats.var_ratio_impr = ostats.var_ratio_impr + delta_ratio * (impr - ostats.avg_ratio_impr);


		}
		ostats.last_quality = quality;

	}
}
/////////////////////////////////////////////////////////////////////////////////////////
// Claspre options
/////////////////////////////////////////////////////////////////////////////////////////
Options::Options()
	: limit(0, 0)
	, timeout(0)
	, fastExit(false) {
}

void Options::initOptions(ProgramOptions::OptionContext& root) {
	using namespace ProgramOptions;
	OptionGroup basic("Basic Options");
	basic.addOptions()
		("time-limit", storeTo(timeout)->arg("<n>"), "Set time limit to %A seconds (0=no limit)")
		("solve-limit", storeTo(limit)->arg("<n>[,<m>]"), "Stop search after <n> conflicts or <m> restarts\n")
		("fast-exit,@1",  flag(fastExit), "Force fast exit (do not call dtors)")
		("file,f,@2", storeTo(input)->composing(), "Input files")
	;
	root.add(basic);
}

bool Options::validateOptions(const ProgramOptions::ParsedOptions& parsed, ProgramOptions::Messages&) {
	if (!parsed.count("solve-limit")) {
		// <TBD>
		// set your default here
		clasp.solve.limit = Clasp::SolveLimits(UINT64_MAX, 4);
	}
	else {
		if (limit.first  == 0) { limit.first = -1; }
		if (limit.second == 0) { limit.second= -1; }
		clasp.solve.limit = Clasp::SolveLimits(static_cast<uint64>(limit.first), static_cast<uint64>(limit.second));
	}
	return true;
}

void Options::applyDefaults(Input::Format f) {
	if (f != Input::SMODELS) {
		Clasp::SatElite::SatElite* pre = new Clasp::SatElite::SatElite();
		pre->options.maxIters = 20;
		pre->options.maxOcc   = 25;
		pre->options.maxTime  = 120;
		clasp.ctx.satPrepro.reset(pre);		
	}
}

bool Options::parsePositional(const std::string&, std::string& out) {
	out = "file";
	return true;
}
/////////////////////////////////////////////////////////////////////////////////////////
// Application
/////////////////////////////////////////////////////////////////////////////////////////
#if !defined(APP_NAME)
#define APP_NAME "claspre"
#endif
#if !defined(CLASPRE_USAGE)
#define CLASPRE_USAGE   APP_NAME " [options] [file]"
#endif
#if !defined (SIGUSR1)
#define SIGUSR1 SIGTERM
#endif
#if !defined(SIGUSR2)
#define SIGUSR2 SIGTERM
#endif

inline bool isStdIn(const std::string& in)  { return in == "-" || in == "stdin"; }
inline bool isStdOut(const std::string& out){ return out == "-" || out == "stdout"; }
/////////////////////////////////////////////////////////////////////////////////////////
// public functions & basic helpers
/////////////////////////////////////////////////////////////////////////////////////////
Application::Application() : timeToFirst_(-1.0), timeToLast_(-1.0), facade_(0), blocked_(0), pending_(0)  {}
Application& Application::instance() {
	static Application inst;
	return inst;
}
void Application::sigHandler(int sig) {
	Application::instance().kill(sig);
}

// Kills any pending alarm
void Application::killAlarm() {
	if (app_.timeout>0) {
		setAlarm(0); 
	}
}

// Called on timeout or signal.
// Prints summary and then kills the application.
void Application::kill(int sig) {
	if (blocked_ == 0) {
		blockSignals();         // ignore further signals
		SCOPE_ALARM_LOCK();
		INFO_OUT(APP_NAME, "INTERRUPTED by signal!");
		if (!facade_ || !facade_->terminate()) {
			if (facade_ && facade_->state() != ClaspFacade::num_states) {
				if (facade_->state() != ClaspFacade::state_start) { timer_[facade_->state()].stop(); }
				timer_[ClaspFacade::state_start].stop();
				cpuTotalTime_.stop();
				printResult(sig);
			}
			appTerminate(app_.clasp.ctx.enumerator()->enumerated > 0 ?  S_SATISFIABLE : S_UNKNOWN);
		}
		else {                  // multiple threads are active - shutdown was initiated
			INFO_OUT(APP_NAME, "Shutting down threads...");
		}
	}
	else if (pending_ == 0) { // signals are currently blocked because output is active
		INFO_OUT(APP_NAME, "Queueing signal...");
		pending_ = sig;
	}
}

// temporarily disable delivery of signals
int Application::blockSignals() {
	return blocked_++;
}

// re-enable signal handling and deliver any pending signal
void Application::unblockSignals(bool deliverPending) {
	if (--blocked_ == 0) {
		int pend = pending_;
		pending_ = 0;
		// directly deliver any pending signal to our sig handler
		if (pend && deliverPending) { kill(pend); }
	}
}

void Application::installSigHandlers() {
	if (signal(SIGINT, &Application::sigHandler) == SIG_IGN) {
		signal(SIGINT, SIG_IGN);
	}
	if (signal(SIGTERM, &Application::sigHandler) == SIG_IGN) {
		signal(SIGTERM, SIG_IGN);
	}
	if (SIGUSR1 != SIGTERM && (signal(SIGUSR1, &Application::sigHandler) == SIG_IGN)) {
		signal(SIGUSR1, SIG_IGN);
	}
	if (SIGUSR2 != SIGTERM && (signal(SIGUSR2, &Application::sigHandler) == SIG_IGN)) {
		signal(SIGUSR2, SIG_IGN);
	}
	if (app_.timeout > 0) {
		setAlarmHandler(&Application::sigHandler);
		if (setAlarm(app_.timeout) == 0) {
			messages.warning.push_back("Could not set time limit!");
		}
	}
}
std::istream& Application::getStream() {
	ProgramOptions::StringSeq& input = app_.input;
	if (input.empty() || isStdIn(input[0])) {
		input.resize(1, "stdin");
		return std::cin;
	}
	else {
		static std::ifstream file;
		if (file.is_open()) return file;
		file.open(input[0].c_str());
		if (!file) { throw std::runtime_error("Can not read from '"+input[0]+"'");  }
		return file;
	}
}

void Application::printHelp(const ProgramOptions::OptionContext& root) {
	printf("%s\n", app_banner);
	printf("usage: %s\n", CLASPRE_USAGE);
	ProgramOptions::FileOut out(stdout);
	root.description(out);
	printf("\n");
	printf(APP_NAME " is part of Potassco: %s\n", "http://potassco.sourceforge.net/");
	printf("Get help/report bugs via : http://sourceforge.net/projects/potassco/support\n");
	fflush(stdout);
}

void Application::printVersion(const ProgramOptions::OptionContext&) {
	printf("%s\n", app_banner);
	printf("\n%s\n%s\n", "Copyright (C) Marius Schneider\n", CLASP_LEGAL);
	fflush(stdout);
}
void Application::printWarnings() const {
	for (ProgramOptions::StringSeq::const_iterator it = messages.warning.begin(); it != messages.warning.end(); ++it) {
		WARNING_OUT(APP_NAME, it->c_str());
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// run - clasp's "main"-function
/////////////////////////////////////////////////////////////////////////////////////////
int Application::run(int argc, char** argv) {
	if (!parse(argc, argv, APP_NAME, &Options::parsePositional)) {
		// command-line error
		ERROR_OUT(APP_NAME, messages.error.c_str());
		INFO_OUT(APP_NAME, "Try '--help' for usage information");
		return S_ERROR;
	}
	if (help || version) {
		return EXIT_SUCCESS;
	}
	installSigHandlers();
	int retStatus = S_UNKNOWN;
	ClaspFacade clasp;
	try {
		StreamInput input(getStream(), detectFormat(getStream()));
		app_.applyDefaults(input.format());
		app_.clasp.ctx.master()->stats.enableStats(2);
		app_.clasp.ctx.master()->stats.enableQueue(100);
		app_.clasp.setMaxSolvers(1);
		app_.clasp.ctx.enableProgressReport(this);
		out_.reset(new Output(app_, clasp));
		facade_ = &clasp;
		cpuTotalTime_.start();
		clasp.solve(input, app_.clasp, this);
		cpuTotalTime_.stop();
		int sig = blockSignals();// disable signal handler
		killAlarm();             // kill any pending alarms;
		printResult(sig);
		if      (clasp.result() == ClaspFacade::result_unsat) retStatus = S_UNSATISFIABLE;
		else if (clasp.result() == ClaspFacade::result_sat)   retStatus = S_SATISFIABLE;
		else                                                  retStatus = S_UNKNOWN;
	}
	catch (const std::bad_alloc&  ) { retStatus = exception(S_MEMORY, "std::bad_alloc"); }
	catch (const std::exception& e) { retStatus = exception(S_ERROR, e.what()); }
	if    (app_.fastExit)           { appTerminate(retStatus);        }
	else                            { fflush(stdout); fflush(stderr); }
	return retStatus;
}

void Application::appTerminate(int status) const {
	fflush(stdout);
	fflush(stderr);
	_exit(status);
}
int Application::exception(int status, const char* what) {
	blockSignals();
	app_.fastExit = true;
	ERROR_OUT(APP_NAME, what);
	if (facade_ && facade_->state() != ClaspFacade::num_states) {
		cpuTotalTime_.stop();
		printResult(status);
	}
	return status;
}
/////////////////////////////////////////////////////////////////////////////////////////
// State & Result functions
/////////////////////////////////////////////////////////////////////////////////////////
// Generates a summary after search has stopped or has been interrupted.
void Application::printResult(int signal) {
	Output::Result r;
	r.signal    = signal;
	r.complete  = (signal == 0 && !facade_->more());
	r.totalTime = timer_[0].total();;
	r.solveTime = timer_[ClaspFacade::state_solve].total();
	r.modelTime = timeToFirst_ != -1.0 ? timeToFirst_ : 0.0;
	double ttl  = timeToLast_ != -1.0 ? timeToLast_ : 0.0;
	r.unsatTime = r.complete && r.solveTime-ttl >= 0.001 ? r.solveTime-ttl : 0.0;
	r.cpuTime   = std::max(cpuTotalTime_.total(), 0.0);
	
	out_->onExit(*app_.clasp.ctx.master(), r);
}

// State-transition callback called by ClaspFacade.
// Handles timing.
void Application::state(ClaspFacade::Event e, ClaspFacade& f) { 
	SCOPE_ALARM_LOCK();
	if (e == ClaspFacade::event_state_enter) {
		timer_[f.state()].start();
	}
	else if (e == ClaspFacade::event_state_exit) {
		timer_[f.state()].stop();
	}
	printWarnings();
	messages.warning.clear();
	out_->onState(e == ClaspFacade::event_state_enter);
}

// Event callback called once by ClaspFacade when the input program
// is ready and once for each model.
void Application::event(const Solver& s, ClaspFacade::Event e, ClaspFacade& f) {
	if (e == ClaspFacade::event_model) {
		timer_[f.state()].lap();
		timeToLast_ = timer_[f.state()].total();
		if (timeToFirst_ == -1.0) {  timeToFirst_ = timeToLast_; }
		const Enumerator& en = *s.sharedContext()->enumerator();
		SCOPE_ALARM_LOCK();
		blockSignals();
		out_->onModel(s, en);
		unblockSignals(true);
	}
	else if (e == ClaspFacade::event_p_prepared) {
		out_->onProgramPrepared(s);
	}
}
void Application::reportProgress(const SolveEvent& e) {
	if (e.type == SolveEvent_t::progress_path) {
		const SolvePathEvent& ev = static_cast<const SolvePathEvent&>(e);
		const Solver& s          = *ev.solver;
		if (ev.evType == SolvePathEvent::event_restart) {
			out_->onRestart(s, ev.cLimit, ev.lLimit);
		}
		else if (ev.type == SolvePathEvent::event_deletion) {
			out_->onDeletion(s, ev.cLimit, ev.lLimit);
		}
	}
}
void Application::reportProgress(const PreprocessEvent& ev) {
	// <TBD>
	// support for preprocessing stats needed?
}


} // end of namespace claspre

