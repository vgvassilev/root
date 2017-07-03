// Author: Enrico Guiraud, Danilo Piparo CERN  03/2017

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RConfigure.h" // R__USE_IMT
#include "ROOT/TDFNodes.hxx"
#include "ROOT/TSpinMutex.hxx"
#include "ROOT/TTreeProcessorMT.hxx"
#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#endif
#include "RtypesCore.h" // Long64_t
#include "TInterpreter.h"
#include "TROOT.h"      // IsImplicitMTEnabled
#include "TTreeReader.h"

#include <cassert>
#include <mutex>
#include <numeric> // std::accumulate
#include <string>
class TDirectory;
class TTree;
using namespace ROOT::Detail::TDF;
using namespace ROOT::Internal::TDF;

namespace ROOT {
namespace Internal {
namespace TDF {

TActionBase::TActionBase(TLoopManager *implPtr, const ColumnNames_t &tmpBranches)
   : fImplPtr(implPtr), fTmpBranches(tmpBranches)
{
}

} // end NS TDF
} // end NS Internal
} // end NS ROOT

TCustomColumnBase::TCustomColumnBase(TLoopManager *implPtr, const ColumnNames_t &tmpBranches, std::string_view name)
   : fImplPtr(implPtr), fTmpBranches(tmpBranches), fName(name){};

ColumnNames_t TCustomColumnBase::GetTmpBranches() const
{
   return fTmpBranches;
}

std::string TCustomColumnBase::GetName() const
{
   return fName;
}

TLoopManager *TCustomColumnBase::GetImplPtr() const
{
   return fImplPtr;
}

TFilterBase::TFilterBase(TLoopManager *implPtr, const ColumnNames_t &tmpBranches, std::string_view name)
   : fImplPtr(implPtr), fTmpBranches(tmpBranches), fName(name){};

TLoopManager *TFilterBase::GetImplPtr() const
{
   return fImplPtr;
}

ColumnNames_t TFilterBase::GetTmpBranches() const
{
   return fTmpBranches;
}

bool TFilterBase::HasName() const
{
   return !fName.empty();
};

void TFilterBase::PrintReport() const
{
   if (fName.empty()) // PrintReport is no-op for unnamed filters
      return;
   const auto accepted = std::accumulate(fAccepted.begin(), fAccepted.end(), 0ULL);
   const auto all = accepted + std::accumulate(fRejected.begin(), fRejected.end(), 0ULL);
   double perc = accepted;
   if (all > 0) perc /= all;
   perc *= 100.;
   Printf("%-10s: pass=%-10lld all=%-10lld -- %8.3f %%", fName.c_str(), accepted, all, perc);
}

// This is an helper class to allow to pick a slot without resorting to a map
// indexed by thread ids.
// WARNING: this class does not work as a regular stack. The size is
// fixed at construction time and no blocking is foreseen.
class TSlotStack {
private:
   unsigned int fCursor;
   std::vector<unsigned int> fBuf;
   ROOT::TSpinMutex fMutex;

public:
   TSlotStack() = delete;
   TSlotStack(unsigned int size) : fCursor(size), fBuf(size) { std::iota(fBuf.begin(), fBuf.end(), 0U); }
   void Push(unsigned int slotNumber);
   unsigned int Pop();
};

void TSlotStack::Push(unsigned int slotNumber)
{
   std::lock_guard<ROOT::TSpinMutex> guard(fMutex);
   fBuf[fCursor++] = slotNumber;
   assert(fCursor <= fBuf.size() && "TSlotStack assumes that at most a fixed number of values can be present in the "
                                    "stack. fCursor is greater than the size of the internal buffer. This violates "
                                    "such assumption.");
}

unsigned int TSlotStack::Pop()
{
   assert(fCursor > 0 &&
          "TSlotStack assumes that a value can be always popped. fCursor is <=0 and this violates such assumption.");
   std::lock_guard<ROOT::TSpinMutex> guard(fMutex);
   return fBuf[--fCursor];
}

TLoopManager::TLoopManager(TTree *tree, const ColumnNames_t &defaultBranches)
   : fTree(std::shared_ptr<TTree>(tree, [](TTree *) {})), fDefaultBranches(defaultBranches),
     fNSlots(TDFInternal::GetNSlots()), fLoopType(ELoopType::kROOTFiles)
{
}

TLoopManager::TLoopManager(Long64_t nEmptyEntries)
   : fNEmptyEntries(nEmptyEntries), fNSlots(TDFInternal::GetNSlots()), fLoopType(ELoopType::kNoFiles)
{
}

/// Run event loop with no source files, in parallel.
void TLoopManager::RunEmptySourceMT()
{
#ifdef R__USE_IMT
   TSlotStack slotStack(fNSlots);
   // Working with an empty tree.
   // Evenly partition the entries according to fNSlots
   const auto nEntriesPerSlot = fNEmptyEntries / fNSlots;
   auto remainder = fNEmptyEntries % fNSlots;
   std::vector<std::pair<Long64_t, Long64_t>> entryRanges;
   Long64_t start = 0;
   while (start < fNEmptyEntries) {
      Long64_t end = start + nEntriesPerSlot;
      if (remainder > 0) {
         ++end;
         --remainder;
      }
      entryRanges.emplace_back(start, end);
      start = end;
   }

   // Each task will generate a subrange of entries
   auto genFunction = [this, &slotStack](const std::pair<Long64_t, Long64_t> &range) {
      auto slot = slotStack.Pop();
      InitAllNodes(nullptr, slot);
      for (auto currEntry = range.first; currEntry < range.second; ++currEntry) {
         RunAndCheckFilters(slot, currEntry);
      }
      slotStack.Push(slot);
   };

   ROOT::TThreadExecutor pool;
   pool.Foreach(genFunction, entryRanges);

#endif // not implemented otherwise
}

/// Run event loop with no source files, in sequence.
void TLoopManager::RunEmptySource()
{
   InitAllNodes(nullptr, 0);
   for (Long64_t currEntry = 0; currEntry < fNEmptyEntries && fNStopsReceived < fNChildren; ++currEntry) {
      RunAndCheckFilters(0, currEntry);
   }
}

/// Run event loop over one or multiple ROOT files, in parallel.
void TLoopManager::RunTreeProcessorMT()
{
#ifdef R__USE_IMT
   TSlotStack slotStack(fNSlots);
   using ttpmt_t = ROOT::TTreeProcessorMT;
   std::unique_ptr<ttpmt_t> tp;
   tp.reset(new ttpmt_t(*fTree));

   tp->Process([this, &slotStack](TTreeReader &r) -> void {
      auto slot = slotStack.Pop();
      InitAllNodes(&r, slot);
      // recursive call to check filters and conditionally execute actions
      while (r.Next()) {
         RunAndCheckFilters(slot, r.GetCurrentEntry());
      }
      slotStack.Push(slot);
   });
#endif // not implemented otherwise
}

/// Run event loop over one or multiple ROOT files, in sequence.
void TLoopManager::RunTreeReader()
{
   TTreeReader r(fTree.get());
   InitAllNodes(&r, 0);

   // recursive call to check filters and conditionally execute actions
   // in the non-MT case processing can be stopped early by ranges, hence the check on fNStopsReceived
   while (r.Next() && fNStopsReceived < fNChildren) {
      RunAndCheckFilters(0, r.GetCurrentEntry());
   }
}

/// Execute actions and make sure named filters are called for each event.
/// Named filters must be called even if the analysis logic would not require it, lest they report confusing results.
void TLoopManager::RunAndCheckFilters(unsigned int slot, Long64_t entry)
{
   for (auto &actionPtr : fBookedActions) actionPtr->Run(slot, entry);
   for (auto &namedFilterPtr : fBookedNamedFilters) namedFilterPtr->CheckFilters(slot, entry);
}

/// Perform clean-up operations. To be called at the end of each event loop.
void TLoopManager::CleanUp()
{
   fHasRunAtLeastOnce = true;

   // forget TActions and detach TResultProxies
   fBookedActions.clear();
   for (auto readiness : fResProxyReadiness) {
      *readiness.get() = true;
   }
   fResProxyReadiness.clear();

   // reset children counts
   fNChildren = 0;
   fNStopsReceived = 0;
   for (auto &ptr : fBookedFilters) ptr->ResetChildrenCount();
   for (auto &ptr : fBookedRanges) ptr->ResetChildrenCount();
   for (auto &pair : fBookedBranches) pair.second->ResetChildrenCount();
}

/// Start the event loop with a different mechanism depending on IMT/no IMT, data source/no data source.
/// Also perform a few setup and clean-up operations (CreateSlots before running, clear booked actions after, etc.).
void TLoopManager::Run()
{
   // jit-compile required actions
   auto error = TInterpreter::EErrorCode::kNoError;
   gInterpreter->ProcessLine(fToJit.c_str(), &error);
   if (error) {
      std::string exceptionText =
         "An error occurred while jitting. The lines above might indicate the cause of the crash\n";
      throw std::runtime_error(exceptionText.c_str());
   }
   fToJit.clear();

   CreateSlots(fNSlots);

#ifdef R__USE_IMT
   if (ROOT::IsImplicitMTEnabled()) {
      switch (fLoopType) {
      case ELoopType::kNoFiles: RunEmptySourceMT(); break;
      case ELoopType::kROOTFiles: RunTreeProcessorMT(); break;
      }
   } else {
#endif // R__USE_IMT
      switch (fLoopType) {
      case ELoopType::kNoFiles: RunEmptySource(); break;
      case ELoopType::kROOTFiles: RunTreeReader(); break;
      }
#ifdef R__USE_IMT
   }
#endif // R__USE_IMT

   CleanUp();
}

/// Build TTreeReaderValues for all nodes
///
/// This method loops over all filters, actions and other booked objects and
/// calls their `BuildReaderValues` methods. It is called once per node per slot, before
/// running the event loop. It also informs each node of the TTreeReader that
/// a particular slot will be using.
void TLoopManager::InitAllNodes(TTreeReader *r, unsigned int slot)
{
   // booked branches must be initialized first
   // because actions and filters might need to point to the values encapsulate
   for (auto &bookedBranch : fBookedBranches) bookedBranch.second->Init(r, slot);
   for (auto &ptr : fBookedActions) ptr->Init(r, slot);
   for (auto &ptr : fBookedFilters) ptr->Init(r, slot);
}

/// Initialize all nodes of the functional graph before running the event loop
///
/// This method loops over all filters, actions and other booked objects and
/// calls their `CreateSlots` methods. It is called once per node before running the
/// event loop. The main effect is to inform all nodes of the number of slots
/// (i.e. workers) that will be used to perform the event loop.
void TLoopManager::CreateSlots(unsigned int nSlots)
{
   for (auto &ptr : fBookedActions) ptr->CreateSlots(nSlots);
   for (auto &ptr : fBookedFilters) ptr->CreateSlots(nSlots);
   for (auto &bookedBranch : fBookedBranches) bookedBranch.second->CreateSlots(nSlots);
}

TLoopManager *TLoopManager::GetImplPtr()
{
   return this;
}

const ColumnNames_t &TLoopManager::GetDefaultBranches() const
{
   return fDefaultBranches;
}

TTree *TLoopManager::GetTree() const
{
   return fTree.get();
}

TCustomColumnBase *TLoopManager::GetBookedBranch(const std::string &name) const
{
   auto it = fBookedBranches.find(name);
   return it == fBookedBranches.end() ? nullptr : it->second.get();
}

TDirectory *TLoopManager::GetDirectory() const
{
   return fDirPtr;
}

void TLoopManager::Book(const ActionBasePtr_t &actionPtr)
{
   fBookedActions.emplace_back(actionPtr);
}

void TLoopManager::Book(const FilterBasePtr_t &filterPtr)
{
   fBookedFilters.emplace_back(filterPtr);
   if (filterPtr->HasName()) {
      fBookedNamedFilters.emplace_back(filterPtr);
   }
}

void TLoopManager::Book(const TmpBranchBasePtr_t &branchPtr)
{
   fBookedBranches[branchPtr->GetName()] = branchPtr;
}

void TLoopManager::Book(const std::shared_ptr<bool> &readinessPtr)
{
   fResProxyReadiness.emplace_back(readinessPtr);
}

void TLoopManager::Book(const RangeBasePtr_t &rangePtr)
{
   fBookedRanges.emplace_back(rangePtr);
}

// dummy call, end of recursive chain of calls
bool TLoopManager::CheckFilters(int, unsigned int)
{
   return true;
}

unsigned int TLoopManager::GetNSlots() const
{
   return fNSlots;
}

/// Call `PrintReport` on all booked filters
void TLoopManager::Report() const
{
   for (const auto &fPtr : fBookedNamedFilters) fPtr->PrintReport();
}

TRangeBase::TRangeBase(TLoopManager *implPtr, const ColumnNames_t &tmpBranches, unsigned int start, unsigned int stop,
                       unsigned int stride)
   : fImplPtr(implPtr), fTmpBranches(tmpBranches), fStart(start), fStop(stop), fStride(stride)
{
}

TLoopManager *TRangeBase::GetImplPtr() const
{
   return fImplPtr;
}

ColumnNames_t TRangeBase::GetTmpBranches() const
{
   return fTmpBranches;
}
