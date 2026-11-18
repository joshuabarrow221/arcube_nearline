////////////////////////////////////////////////////////////////////////////////////////////
// ORIGINAL SCRIPT DEVELOPED BY **GIANFRANCO INGRATTA** (York University, ingratta@yorku.ca)
// WITH SUPPORT FROM SUPERVISOR BRUCE HOWARD (York University, blhoward@yorku.ca)
////////////////////////////////////////////////////////////////////////////////////////////

// Script updated by J. L. Barrow (UMN, jbarrow@umn.edu)
// with **much help** from Anthropic Claude Fable 5 High effort.
//
// Reconciled line-by-line against NOvA's current IFDBSpillInfo module and
// CAFAna spill cuts (novasoft main @ d4b6064706, 2026-11-12) by
// J. L. Barrow with help from Anthropic Claude Opus 5.5 (Nov. 2026). The
// same algorithms are implemented for CAF production in DUNE's
// ND_CAFMaker (src/beam/NuMIBeamQuality.{h,cxx}); the two are meant to
// agree spill by spill and were cross-checked against each other on
// 2x2 Run 1 data.

// =====================================================================
//  get_data.cpp
//
//  Fetch per-spill NuMI beam data for a given time window from the
//  IFBeam database and write out
//    (1) a ROOT TTree "beam" + a tidy long-format CSV: every spill,
//        deliberately UNCUT, every device reading as archived;
//    (2) a parallel TTree "quality" + *_quality.csv: NOvA's per-spill
//        "goodbeam" evaluation (derived quantities, the individual
//        pass flags and the combined verdict);
//    (3) a record of the exact cut values used (TParameter/TNamed
//        objects in the ROOT file and a *_cuts.json sidecar), which the
//        plotting/report scripts in this directory read so that every
//        figure draws the cuts that were actually applied.
//
//  Intended use: associating DUNE ND 2x2 demonstrator triggers with the
//  NuMI beam conditions of the corresponding spill (protons on target,
//  horn current and polarity, beam position and width at the target).
//
//  ---------------------------------------------------------------
//  ACNET vs. the IFBeam database — which is which?
//  ---------------------------------------------------------------
//  The device names used below ("E:TRTGTD", "E:NSLINA", ...) are ACNET
//  device identifiers. ACNET is Fermilab's accelerator control network;
//  every instrument in the accelerator complex (toroids, BPMs, power
//  supplies, ...) is addressed by such a name, with the leading letter
//  ("E:") denoting the ACNET device class/area (external-beamline
//  devices, which includes the NuMI primary proton line).
//
//  This program, however, never talks to ACNET directly. It queries the
//  IFBeam database (IFDB, "Intensity Frontier Beam Database"): an
//  archive maintained for the IF experiments (MINOS, MINERvA, NOvA,
//  DUNE/2x2, ...) into which a curated subset of ACNET device readings
//  is copied, timestamped by accelerator timing (TCLK) events. IFDB is
//  exposed through a public, unauthenticated HTTPS REST API — no
//  Kerberos, VPN, or Fermilab credentials are needed. So: the
//  *variables* are ACNET devices; the *source we read them from* is the
//  IFBeam database's archive of those devices. (NOvA reads the same
//  archive through the "NuMI_Physics_A9" bundle, which holds 65 devices
//  per $A9 record; the 15 fetched here are the ones its goodbeam
//  selection needs, plus the TOR101 direct readback and E:HRNDIR.)
//
//  ---------------------------------------------------------------
//  Accelerator clock events: TCLK $A9 (archive) vs. MIBS $74 (beam sync)
//  ---------------------------------------------------------------
//  Each fetched reading is tagged with the timing event on which it was
//  recorded. NOvA's spill-server message types (novasoft
//  Commissioning/Nearline/TimingPlots/SpillMessage.h) distinguish:
//      kNuMI            "MIBS $74 proton extraction into NuMI"
//      kNuMItclk        "tevatron clock, either $A9 or $AD depending on
//                        xml parameter"
//      kNuMISampleTrig  "$A4, NuMI cycle sample trigger, reference for $A5"
//      kNuMIReset       "$A5, NuMI reset for beam"
//      kSuperCycle      "$00, Super cycle and master clock reset"
//  NOvA's detectors are triggered from the MIBS $74 beam-sync event,
//  whereas the IFBeam archive files every NuMI device reading under the
//  TCLK $A9 event, with a millisecond "clock" timestamp. We therefore
//  request event "$A9" (query string "e=e,a9") and obtain exactly one
//  reading per NuMI spill rather than the continuous datalogger stream.
//  A direct check (2026-11-18): E:TOR101 over 2024-07-09 00:54-00:55
//  CDT returns 46 records on "e,a9" and none on "e,74", "e,ad" or
//  "e,a5" — the archive holds no separate $74 record. A detector
//  trigger is therefore associated with its $A9 record (NOvA requires
//  |t_trigger - t_IFDB| < 0.5 s, cut 6 below), and the trigger-minus-$A9
//  offset measured in data is itself evidence of which accelerator event
//  drives a detector's beam trigger. In July 2024 the $A9 records arrive
//  every 1.200 s, with one ~6.47 s gap per ~60.46 s Main Injector
//  super-cycle (the slot given to the Switchyard slow extraction; not
//  data loss — no NuMI spill exists in that slot).
//
//  ---------------------------------------------------------------
//  References (publicly accessible unless noted)
//  ---------------------------------------------------------------
//  [1] P. Adamson et al., "The NuMI Neutrino Beam",
//      Nucl. Instrum. Meth. A 806 (2016) 279-306, arXiv:1507.06690,
//      https://arxiv.org/abs/1507.06690
//      Sec. 4 describes the primary-beam instrumentation quoted below:
//      the two POT toroids (TOR101 upstream in the NuMI line, TORTGT
//      ~10 m before the target), the capacitive beam position monitors
//      (BPMs), the profile-monitor multiwires, and the focusing horns
//      with forward (FHC, neutrino) / reverse (RHC, antineutrino)
//      current operation.
//  [2] NOvA, IFDBSpillInfo/IFDBSpillInfo_module.cc and IFDBSpillInfo.fcl
//      (novasoft main) — the de-facto reference implementation for
//      turning these IFDB devices into per-spill beam-quality
//      quantities. Public mirror:
//      https://github.com/novaexperiment/novasoft/blob/main/IFDBSpillInfo/IFDBSpillInfo_module.cc
//      (historical Redmine link, now behind Fermilab SSO:
//      https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/IFDBSpillInfo_module.cc#L685).
//      Credit where due: the module was originally authored by
//      L. Goodenough (ANL) and S. Phan-Budd (Winona State U.); the
//      intensity-weighted extrapolated-to-target position machinery
//      (BpmProjection/BpmAtTarget) ported below is due to G. Brunetti
//      (FNAL); T. Lackey (FNAL) is its long-term maintainer
//      (period-dependent cuts, robust missing-toroid handling);
//      R. Schroeter added the is0HC horn-off flag; A. Radovic widened
//      the position cuts for new beam modes and introduced the
//      E:HRNDIR-based RHC flag; J. Paley the RunHistory RHC
//      determination; the stripline calibration constants are from
//      J. Hylen (FNAL AD).
//  [3] NOvA, CAFAna/Cuts/SpillCuts.cxx (novasoft main): the analysis-
//      level copy of the same cuts (kTightBeamQualityCuts,
//      kBeamPositionCut, kBeamWidthCut, kHornPolarityCheck), which
//      carries the NOvA-docdb provenance restated below.
//  [4] IFBeam REST API syntax, "DataAccessSyntax" wiki of the
//      "ifbeamdata" project,
//      https://cdcvs.fnal.gov/redmine/projects/ifbeamdata/wiki/DataAccessSyntax
//      (now behind Fermilab SSO; snapshots exist on web.archive.org).
//      Documents the query parameters used in createUrl(): v=device,
//      e=event, t0/t1=time window, f=output format.
//
//  ---------------------------------------------------------------
//  Beam-quality ("goodbeam") cuts: NOvA's criteria and this code
//  ---------------------------------------------------------------
//  NOvA's per-spill beam-quality selection is the goodbeam block at the
//  end of IFDBSpillInfo::produce [2], with numerical values from
//  IFDBSpillInfo.fcl (standard_ifdbspillinfo) and an identical copy in
//  CAFAna's kTightBeamQualityCuts [3]. A spill is "goodbeam" only if
//  ALL of the following hold (strict inequalities, as in [2]; CAFAna's
//  copy uses inclusive bounds, which differ only at exact equality):
//
//    1. POT:          spillpot > 2.00e12 protons (MinPOTCut, "x10^12
//                     pot"). POT is E:TRTGTD, falling back to E:TR101D
//                     when TRTGTD is missing or reads < 0.02e12;
//                     negative toroid readings are clamped to 0 ("This
//                     guarantees that when toroid values go negative,
//                     total POT is not reduced").
//    2. Horn current: -202 kA < hornI < -196.4 kA (Min/MaxHornICut),
//                     where hornI is the calibrated sum of E:NSLINA-D
//                     (see the NSLIN block below). Spills with
//                     |hornI| < 1 kA (HornICut0HC) are flagged horn-off
//                     (is0HC). This window applies to FHC AND RHC
//                     running alike: see "Horn polarity" below.
//    3. Position x:   -2.00 mm < posx < +2.00 mm (Min/MaxPosXCut).
//    4. Position y:   epoch dependent (Min/MaxPosYCut vectors, table
//                     below); (-3.50, +0.50) mm for all running since
//                     NOvA period 13.
//                     NOTE: NOvA's posx/posy are NOT the raw E:H/VPTGT
//                     readbacks: they are the intensity-weighted
//                     per-batch positions (BPM array element 0, the
//                     auto-tune average, skipped; batches with zero BPM
//                     intensity skipped) linearly extrapolated from the
//                     121 and TGT BPM stations to the target z
//                     (BpmProjection/BpmAtTarget in [2]).
//    5. Width:        0.57 mm < widthx, widthy < max (epoch dependent:
//                     1.58 or 1.88 mm, table below), where the widths
//                     are Gaussian sigmas fitted to the E:MTGTDS[]
//                     multiwire profile (horizontal wires at array
//                     indices 103-150, vertical at 151-198;
//                     ProfileProjection/GetGaussFit in [2]).
//    6. Timing:       |trigger time - IFDB spill time| < 0.5 s
//                     (MaxDeltaTCut = 0.5e9 ns): "This checks that the
//                     event time is close to the database time to
//                     ensure that we are grabbing the right data" [2].
//
//  NOvA selects the position-y and width cuts by detector run number
//  ("Periods 1-10 / 11-12 / 13-?" in [2] and [3]):
//
//     NOvA period  ND runs       FD runs       posy window [mm]  max widthx  max widthy
//     1-9          < 13428       < 34574       (-2.00, +2.00)    1.58        1.58
//     10           13428-13886   34574-38300   (-2.00, +2.00)    1.88        1.58
//     11-12        13887-14831   38301-44836   (-4.50,  0.00)    1.88        1.58
//     13+          >= 14832      >= 44837      (-3.50, +0.50)    1.88        1.88
//
//  Period 13 began in late 2022, so EVERY 2x2 data-taking period
//  (Run 1 = July 2024 onward) falls in the "13+" row, which is this
//  code's default (novaEpoch = "P13+"). Earlier epochs remain selectable
//  for completeness. (An earlier version of this file applied the
//  periods 1-10 values — posy within +-2 mm, widths below 1.58 mm — and
//  described the order of the epochs backwards. On the July 2024 Run 1
//  windows the change alters 0 of 8,443 spill verdicts, because the
//  beam sat at posy -1.65 to -0.33 mm with widths 1.42-1.54 mm.)
//
//  ---------------------------------------------------------------
//  Where NOvA's cut values come from (restated from novasoft)
//  ---------------------------------------------------------------
//  Restated from the comments in [2] and [3] and from the commit
//  messages that introduced each value (git log of IFDBSpillInfo/ on
//  novasoft main):
//
//   2014-06-05  L. Goodenough    IFDBSpillInfo added, "which will replace
//                                the antiquated NuMISpillInfo Package".
//   2014-12-01  L. Goodenough    "made changes to beam parameter cuts based
//                                on discussions in Beam Group. changes to
//                                horn current cuts are relatively large,
//                                but the range of horn currents measured
//                                thus far are well within the new range
//                                of -202kA to -198kA. the delta t cut
//                                ensures that we are matching the spill
//                                with the correct IFDB data."
//   2015-05-16  G. Davies        IFBeam cache TimeWindow 120 s -> 600 s.
//   2015-11-06  A. Radovic       "Increased range of MinPosYCut and
//                                MinPosXCut to allow new beam mode running
//                                to pass. We will want a more robust
//                                version of this ASAP."
//   2016-09-02  R. Schroeter     "support for new is0HC (horn off)
//                                spillinfo variable."
//   2020-03-10  R. Bowles        "Changed the maximum X Beam Width from
//                                1.58mm to 1.88mm" (period 10 onward;
//                                [3]: "For a summary of the reasoning,
//                                see ... docid=44208", NOvA docdb-44208).
//   2021-01-12  T. Lackey        "Beam position moved down 2 mm for period
//                                11. Changes will use old cuts for any
//                                FD/ND runs prior to period 11."
//   2021-01-27  J. M. Carceller  "Change cut on the min Y position of the
//                                beam from -4 to -4.5" (periods 11-12).
//   2022-01-10  T. Lackey        "Updated horn current max cut from -198
//                                to -196.4 kA" (novasoft PR #166); [3]:
//                                "Max cut changed from -198. Value from
//                                Katsuya. More details in docdb 49148
//                                slide 15".
//   2022-11-22  T. Lackey        "Updated vertical beam position cuts for
//                                period 13. More details in docdb-57045."
//   2022-12-01  T. Lackey        "Update max y-width cut to 1.88 mm for
//                                period 13"; [3]: "Apply different y
//                                width cuts starting in period 13. See
//                                docdb-57045, talk from 2022-12-05."
//   2022-12-21  T. Lackey        "Implement period-dependant max beam
//                                width cuts in NearlineAna and
//                                IFDBSpillInfo modules."
//   2024-03-29  T. Lackey        Split the single IFBeam try/catch "for
//                                events where E:TRTGTD is fine but
//                                E:TR101D is missing" so the other
//                                quantities are still filled; "Remove
//                                horndir since it's not used anywhere; we
//                                use RunHistory to determine RHC status
//                                instead."
//   2025-08-22  T. Lackey        gpsspilltime{sec,nsec} now carry the
//                                NOvA DAQ (detector) event time.
//
//  [3] also keeps an older, looser selection (kLooseBeamQualityCuts,
//  "docdb 12474 slides 3-10": |dt| <= 1 s, POT >= 0.5e12,
//  -205 <= hornI <= -195 kA, same position/width cuts). The horn
//  stripline calibration is attributed in [2] to "Updated numbers from
//  Jim Hylen, email Oct 4, 2013".
//
//  ---------------------------------------------------------------
//  Horn polarity: E:HRNDIR, not the sign of the horn current
//  ---------------------------------------------------------------
//  The calibrated stripline sum is ~ -200 kA in BOTH forward (FHC,
//  neutrino-mode) and reverse (RHC, antineutrino-mode) running: NOvA
//  applies the same -202 < hornI < -196.4 kA window to all spills and
//  takes the polarity from its run database (RunHistory::IsRHC) [2].
//  The history in [2]: on 2016-06-28 A. Radovic briefly cut on
//  |hornI| "to prevent the good spill cut from breaking during RHC
//  running", then on 2016-06-30 restored the signed window and "Uses the
//  horndir beam data to fill the isRHC spillinfo" (isRHC = horndir > 0);
//  on 2018-04-16 J. Paley moved isRHC to RunHistory run ranges, "instead
//  of the unreliable horn current readings", and the horndir code was
//  removed in 2024. Checked directly against the IFBeam archive
//  (2026-11-18, 1-minute $A9 samples at the dates below; epochs from the
//  RHC run ranges documented in [3]):
//
//     date (CDT)         NOvA epoch                 E:HRNDIR    E:NSLINA
//     2017-04-12 12:00   RHC (range 3, epoch 6a)    +742 V      -49.6 kA
//     2018-03-15 12:00   RHC (range 4, epochs 7)    +742 V      -49.5 kA
//     2019-02-15 12:00   RHC (range 5, epochs 8)    +742 V      -49.5 kA
//     2021-06-15 12:00   FHC                        -748 V      -49.1 kA
//     2023-03-15 12:00   FHC                        -753 V      -49.4 kA
//     2024-07-09 00:54   RHC (range 6; 2x2 Run 1)   +19.4 V     -49.7 kA
//
//  So E:HRNDIR > 0 <=> RHC throughout, and the stripline (hence hornI)
//  sign does NOT change with polarity. Note the much smaller RHC readback
//  since the 2023-24 shutdown (+19 V instead of ~+742 V): the sign is
//  still unambiguous, but the magnitude should not be relied upon, and
//  the convention should be re-checked whenever a new running period
//  starts. (A horn-off sample, 2015-06-15, read NSLINA = 0 and
//  E:HRNDIR = -1.4 V: polarity is meaningless when is0HC.) [3] (novasoft
//  main, commit 0acaf7b057 of 2026-09-16) lists NOvA's sixth RHC range as
//  "from long post-shutdown 2023-24. (epoch 15a)": ND runs 15753-15922 and
//  FD runs 49841-50993 (RunHistory.fcl RHC stop runs 15922 and 50993). The
//  2x2 Run 1 data (July 2024, E:HRNDIR = +19 V) fall inside it and are
//  antineutrino-mode (RHC) data. NOvA's run lists mark no later running
//  as RHC (i.e. FHC or horn-off afterwards).
//  No horn-on IFBeam sample after it was found to re-check the sign
//  convention (2024-12 to 2026-05: no records in the sampled minutes;
//  2026-10-07: horn off, NSLINA = 0, E:HRNDIR = -1.9 V), so check it again
//  on the first horn-on data of a new running period.
//  (An earlier version of this file called the horn-current window an
//  "FHC selection"; it is not — it is NOvA's selection for both
//  polarities.)
//
//  ---------------------------------------------------------------
//  What THIS code applies
//  ---------------------------------------------------------------
//  The "beam" tree and main CSV remain deliberately UNCUT — every spill
//  returned by IFBeam is written out so downstream studies can see
//  off-nominal beam. The parallel "quality" tree EVALUATES cuts 1-5 per
//  spill exactly as NOvA does — BpmProjection / BpmAtTarget /
//  ProfileProjection / GetStats / GetGaussFit of [2] are ported below,
//  including one NOvA quirk kept on purpose for spill-by-spill
//  consistency (see bpmAtTarget()) — and stores the derived quantities,
//  the six pass flags and the combined goodbeam verdict. Nothing is
//  filtered: failing spills are written with goodbeam = 0. Cut 6
//  (trigger timing) is NOT evaluated — it compares a detector trigger
//  time against the spill record and this standalone fetcher has no
//  detector trigger stream — so "goodbeam" here means cuts 1-5 (the
//  ND_CAFMaker implementation evaluates all six). The only other
//  selections made here are data-hygiene:
//    - readings from different devices are associated with the same
//      spill only if their timestamps agree within beamMatchDT = 0.5 s
//      (the same half-width as NOvA's timing cut 6);
//    - the derived E:NSLIN horn current is computed only for spills
//      where all four striplines have a time-matched record (NOvA reads
//      the four in a single call and gets nothing if one is missing);
//    - unmatched devices appear as NaN (scalars) or empty vectors
//      (arrays) in the TTree — nothing is silently dropped.
//  Where NOvA would store a sentinel (-9 m, i.e. -9000 mm, or 0) for an
//  undefined quantity, this code stores NaN; the pass/fail verdicts are
//  identical because both fail every window.
//
//  The plotting scripts (beam_gif.py, beam_batches_gif.py,
//  beam_slosh.C) plot each batch's 121- and TGT-station readings
//  linearly extrapolated to the target z with the same geometry as
//  bpmProjection() below and DRAW the position window of cuts 3-4 that
//  is recorded in *_cuts.json / the ROOT file, for illustration only;
//  they differ from the cut quantity only in skipping the intensity
//  weighting (individual batches, not the per-spill weighted mean).
// =====================================================================

#include <curl/curl.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <TNamed.h>
#include <TParameter.h>

using json = nlohmann::json;

// One device's archive for the queried window:
//   key   = spill time, seconds since the Unix epoch (IFBeam "clock"),
//   value = the reading(s) on that spill — a single number for scalar
//           devices (toroids, horn striplines), an array for the
//           BPM ("[]") and multiwire devices.
using BeamInfo = std::map<double, std::vector<double>>;
using DeviceMap = std::map<std::string, BeamInfo>;

// Half-width (seconds) of the coincidence window used to decide that
// readings from two different devices belong to the same spill. Spills
// are >~1.2 s apart, so +-0.5 s matches at most one spill.
double beamMatchDT = 0.5;

// =====================================================================
// NOvA goodbeam cut values, per NOvA running epoch.
// =====================================================================
// All numbers are IFDBSpillInfo.fcl's standard_ifdbspillinfo values [2],
// converted to this file's units (absolute protons, kA, mm, s). The
// provenance of every value is restated in the header comment above.
struct NOvABeamCuts {
    std::string epoch;            // NOvA period label these values belong to
    double minPOT;                // protons; goodbeam needs spillpot > minPOT
    double minHornI, maxHornI;    // kA; needs minHornI < hornI < maxHornI
    double hornI0HC;              // kA; |hornI| < hornI0HC flags horn-off (0HC)
    double minPosX, maxPosX;      // mm at target
    double minPosY, maxPosY;      // mm at target
    double minWidthX, maxWidthX;  // mm (multiwire Gaussian sigma)
    double minWidthY, maxWidthY;  // mm
    double maxDeltaT;             // s; |t_trigger - t_IFDB| < maxDeltaT (cut 6)
};

// Cut values for one NOvA epoch label:
//   "P1-9"   NOvA periods 1-9    (ND runs < 13428, FD runs < 34574)
//   "P10"    NOvA period 10      (ND 13428-13886, FD 34574-38300)
//   "P11-12" NOvA periods 11-12  (ND 13887-14831, FD 38301-44836)
//   "P13+"   NOvA periods 13+    (ND >= 14832, FD >= 44837) — all 2x2 data
// These mirror the run-number branches in IFDBSpillInfo::produce [2] and
// in CAFAna's kBeamPositionCut / kBeamWidthCut [3].
NOvABeamCuts novaCutsForEpoch(const std::string& epoch) {
    NOvABeamCuts c;
    c.epoch = epoch;

    // ---- epoch-independent values ----
    c.minPOT    = 2.00e12;  // MinPOTCut: 2.00 "x10^12 pot" -> absolute protons
    c.minHornI  = -202.0;   // MinHornICut [kA]
    c.maxHornI  = -196.4;   // MaxHornICut [kA] (-198 before 2022-01-10, docdb 49148)
    c.hornI0HC  = 1.0;      // HornICut0HC [kA]
    c.minPosX   = -2.00;    // MinPosXCut [mm]
    c.maxPosX   = +2.00;    // MaxPosXCut [mm]
    c.minWidthX = 0.57;     // MinWidthXCut [mm]
    c.minWidthY = 0.57;     // MinWidthYCut [mm]
    c.maxDeltaT = 0.5;      // MaxDeltaTCut: 0.5e9 ns -> s

    // ---- epoch-dependent values: MinPosYCut/MaxPosYCut and
    //      MaxWidthXCut/MaxWidthYCut are vectors in the fcl, indexed by
    //      period exactly as below ----
    if (epoch == "P1-9") {
        c.minPosY = -2.00; c.maxPosY = +2.00;      // index 0
        c.maxWidthX = 1.58; c.maxWidthY = 1.58;    // index 0, index 0
    } else if (epoch == "P10") {
        c.minPosY = -2.00; c.maxPosY = +2.00;      // index 0
        c.maxWidthX = 1.88; c.maxWidthY = 1.58;    // index 1, index 0
    } else if (epoch == "P11-12") {
        c.minPosY = -4.50; c.maxPosY = 0.00;       // index 1 ("moved down 2 mm")
        c.maxWidthX = 1.88; c.maxWidthY = 1.58;    // index 1, index 0
    } else if (epoch == "P13+") {
        c.minPosY = -3.50; c.maxPosY = +0.50;      // index 2 (docdb-57045)
        c.maxWidthX = 1.88; c.maxWidthY = 1.88;    // index 1, index 1 (docdb-57045)
    } else {
        throw std::invalid_argument("novaCutsForEpoch: unknown NOvA epoch '" + epoch +
                                    "' (use P1-9, P10, P11-12 or P13+)");
    }
    return c;
}

// ---------------------------------------------------------------------
// URL construction for the IFBeam REST API [4].
//
// Endpoint: https://dbdata3vm.fnal.gov:9443/ifbeam/data/data
// Query parameters:
//   v  = ACNET device name; a trailing "[]" requests the full array
//        readback rather than only the first element,
//   e  = timing event on which the reading was archived; "e,a9" is the
//        TCLK NuMI event $A9 (one reading per spill; see header),
//   t0 = window start, ISO-8601 with UTC offset,
//   t1 = window end,
//   f  = response format ("json" here; csv/xml also exist).
// ---------------------------------------------------------------------
std::string createUrl(const std::string& device,
                      const std::string& min_time_iso,
                      const std::string& max_time_iso) {
    // Device names (e.g. "E:HPTGT[]") and ISO times ("+", ":") contain
    // characters that must be percent-encoded in a query string.
    CURL* curl = curl_easy_init();
    if (!curl) {
        throw std::runtime_error("createUrl: failed to initialize curl for escaping");
    }
    auto escape = [curl](const std::string& s) {
        char* e = curl_easy_escape(curl, s.c_str(), static_cast<int>(s.size()));
        std::string out = e ? e : s;
        curl_free(e);
        return out;
    };
    std::ostringstream url_stream;
    url_stream << "https://dbdata3vm.fnal.gov:9443/ifbeam/data/data?v=" << escape(device)
               << "&e=" << "e,a9"
               << "&t0=" << escape(min_time_iso)
               << "&t1=" << escape(max_time_iso)
               << "&f=json";
    curl_easy_cleanup(curl);
    return url_stream.str();
}

// ---------------------------------------------------------------------
// Unit-string -> scale factor.
//
// IFBeam stores each reading's engineering unit as a string, and for
// the toroids that string *is* the power of ten, e.g. "E12": a TRTGTD
// value of 49.77 with units "E12" means 4.977e13 protons. This helper
// turns "E<n>" into 10^n and leaves genuine units (e.g. "mm") at 1.0.
// (NOvA reads the raw values through ifbeam's BeamFolder and therefore
// keeps POT "normalized by 10^12"; here POT is absolute everywhere.)
// ---------------------------------------------------------------------
double unitToFactor(const std::string& unit) { //because the exponent part in IFBeam is stored as a string under "value" :/
      std::regex re("E(\\d+)");
      std::smatch match;

      if (std::regex_search(unit, match, re) && match.size() > 1) {
          int exponent = std::stoi(match.str(1));
          return std::pow(10, exponent);
      } else {
          return 1.0;
      }
  }

// ---------------------------------------------------------------------
// Epoch (seconds, possibly with a fractional part) -> ISO8601 conversion,
// in the local time zone of the machine, with an explicit UTC offset.
// ---------------------------------------------------------------------
std::string toISO8601(double time_sec) {
    time_t seconds = static_cast<time_t>(std::floor(time_sec));
    long long ms = static_cast<long long>(std::llround((time_sec - static_cast<double>(seconds)) * 1000.0));
    if (ms >= 1000) { seconds += 1; ms -= 1000; }   // rounding carried into the next second

    // std::localtime is NOT thread-safe (it uses a shared static buffer).
    // We call the re-entrant variant once and reuse the result.
    std::tm local_tm{};
#if defined(_WIN32)
    localtime_s(&local_tm, &seconds);
#else
    localtime_r(&seconds, &local_tm);
#endif

    std::ostringstream ss;
    ss << std::put_time(&local_tm, "%Y-%m-%dT%H:%M:%S");
    if (ms > 0) {
        ss << '.' << std::setfill('0') << std::setw(3) << ms;
    }

    long gmtoff = local_tm.tm_gmtoff;
    ss << (gmtoff >= 0 ? '+' : '-')
       << std::setfill('0') << std::setw(2) << std::abs(gmtoff) / 3600
       << ':' << std::setfill('0') << std::setw(2) << (std::abs(gmtoff) % 3600) / 60;

    return ss.str();
}

// ---------------------------------------------------------------------
// curl write callback: append the response body to a std::string.
// ---------------------------------------------------------------------
static size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* userp) {
    userp->append(static_cast<char*>(contents), size * nmemb);
    return size * nmemb;
}

// ---------------------------------------------------------------------
// Query one device from the IFBeam database and unpack the JSON.
//
// The response is a Google-Charts-style table {"cols": [...], "rows":
// [...]}, one row per spill, each row holding:
//   "clock" — the reading's timestamp in *milliseconds* since the Unix
//             epoch (converted to seconds here),
//   "units" — engineering-unit string, see unitToFactor(),
//   "c"     — the cell array: c[0] is the time cell (a "Date(...)"
//             string, skipped because it is not a number), c[1..] are
//             the readings (1 for scalar devices; 7 for the BPM arrays;
//             216 for the MTGTDS multiwire).
// ---------------------------------------------------------------------
BeamInfo retrieveInfoFromDataBase(const std::string& url) {
    constexpr double ms_to_s = 1e-3;
    BeamInfo data;

    CURL* curl = curl_easy_init();
    if (!curl) {
        throw std::runtime_error("retrieveInfoFromDataBase: curl_easy_init() failed");
    }

    std::string readBuffer;
    long httpCode = 0;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 60L);

    CURLcode res = curl_easy_perform(curl);
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &httpCode);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK) {
        // Throw rather than abort: the caller decides whether to retry,
        // skip the device, or terminate.
        throw std::runtime_error("curl_easy_perform() failed: " +
                                 std::string(curl_easy_strerror(res)));
    }
    if (httpCode != 200) {
        // e.g. HTTP 400 for a malformed query or 500 for a server error;
        // the body is then an HTML error page, not JSON.
        throw std::runtime_error("IFBeam returned HTTP " + std::to_string(httpCode) +
                                 " for url '" + url + "'");
    }

    try {
        auto json_data = json::parse(readBuffer);
        for (const auto& spill : json_data["rows"]) {
            double time = spill["clock"].get<double>() * ms_to_s;
            std::string unit = spill.value("units", "");
            std::vector<double> values;

            if (spill.contains("c")) {
                for (const auto& c_entry : spill["c"]) {
                    if (c_entry.contains("v") && c_entry["v"].is_number()) {
                        values.push_back(c_entry["v"].get<double>() * unitToFactor(unit));
                    }
                }
            }
            data[time] = std::move(values);
        }
    } catch (const json::exception& e) {
        throw std::runtime_error("JSON parsing error for url '" + url + "': " + e.what());
    }

    return data;
}

// ---------------------------------------------------------------------
// Find the record of `info` whose spill time is within `dt` seconds of
// `t`, as an iterator (end() if none). Spills are >= 1.2 s apart, so
// with dt = 0.5 s at most one record can lie in the window.
// ---------------------------------------------------------------------
BeamInfo::const_iterator findMatchIt(const BeamInfo& info, double t, double dt) {
    auto it = info.lower_bound(t - dt);
    if (it != info.end() && std::fabs(it->first - t) <= dt) return it;
    return info.end();
}

// Same, returning a pointer to the matched reading(s) or nullptr.
const std::vector<double>* findMatch(const BeamInfo& info, double t, double dt) {
    auto it = findMatchIt(info, t, dt);
    return (it != info.end()) ? &it->second : nullptr;
}

// =====================================================================
// NOvA-style beam-quality computation, ported from IFDBSpillInfo [2].
// =====================================================================

// BPM station z-positions [feet], from [2]: "new beam monitor positions
// after survey. Origin is between two budal monitors of the NOVA target"
// (pre-survey values, kept there as comments: HP121 -72.2309, VP121
// -71.3142, HPTGT -33.1564, VPTGT -32.2397 ft). Only ratios enter the
// linear extrapolation, so the unit is irrelevant; NOvA's z_targ = 0
// ("In NOvA era (for now) we will assume the following. In MINOS era,
// target location was 3937 after time period 1155564632").
constexpr double kZhp121 = -68.04458, kZvp121 = -66.99283,
                 kZhptgt = -31.25508, kZvptgt = -30.16533, kZtarg = 0.0;
// "wire spacing on SWIC" [2]
constexpr double kWireSpacing = 0.5;   // mm

// Port of extrapolate_position [2]: "Linear extrapolation to location
// z3" of a transverse position t measured at z1 and z2.
double extrapolatePosition(double t1, double z1, double t2, double z2, double z3) {
    return t1 + (t2 - t1) * (z3 - z1) / (z2 - z1);
}

// Port of BpmProjection [2]. Inputs are the six BPM arrays of one spill
// (7 elements each). For every Booster batch it returns
//   xp, yp : the horizontal/vertical beam position linearly extrapolated
//            from the 121 to the TGT station and on to the target z [mm],
//   xi, yi : the HITGT/VITGT intensities ("Simply transfer intensities
//            over").
// As in NOvA: the array length used is the shortest of the six; element
// 0 is skipped "(average used in auto-tune) but only if there are
// multiple batches"; and for batches where either target-BPM intensity
// reads exactly 0 the position AND intensity entries are set to 0
// ("If they are, set position values to 0"), which makes BpmAtTarget
// skip them. Positions are in mm here (NOvA works in metres and
// multiplies by 1000 when filling bposx/bposy).
void bpmProjection(const std::vector<double>& vp121, const std::vector<double>& hp121,
                   const std::vector<double>& vptgt, const std::vector<double>& hptgt,
                   const std::vector<double>& vitgt, const std::vector<double>& hitgt,
                   std::vector<double>& xp, std::vector<double>& yp,
                   std::vector<double>& xi, std::vector<double>& yi) {
    xp.clear(); yp.clear(); xi.clear(); yi.clear();
    // "Return shortest array length of all BPM positions and intensities used"
    const size_t n = std::min({vp121.size(), hp121.size(), vptgt.size(),
                               hptgt.size(), vitgt.size(), hitgt.size()});
    const size_t start = (n > 1) ? 1 : 0;   // skip auto-tune element 0
    for (size_t ind = start; ind < n; ++ind) {
        // "Check to see if elements of BPMS[4] and BPMS[5] (BPM intensities) are 0."
        if (vitgt[ind] == 0.0 || hitgt[ind] == 0.0) {
            yp.push_back(0.0); xp.push_back(0.0);
            yi.push_back(0.0); xi.push_back(0.0);
            continue;
        }
        // "Returns transverse position (x or y) at z_targ assuming a
        //  linear extrapolation from z_vp121 to z_vptgt"
        yp.push_back(extrapolatePosition(vp121[ind], kZvp121, vptgt[ind], kZvptgt, kZtarg));
        xp.push_back(extrapolatePosition(hp121[ind], kZhp121, hptgt[ind], kZhptgt, kZtarg));
        // "Returns the horizontal or vertical intensity of beam"
        yi.push_back(vitgt[ind]);
        xi.push_back(hitgt[ind]);
    }
}

// Port of BpmAtTarget [2]: intensity-weighted mean of the per-batch
// target positions from bpmProjection(). Batches with a zero intensity
// in either plane are skipped ("If average intensity is 0, skip to next
// value") and negative intensities count as zero ("if I>0, set I=I,
// else set I=0"). Returns false — leaving posx/posy untouched — when no
// positive total intensity remains (NOvA then keeps its -9 m default).
//
// NOvA quirk, KEPT ON PURPOSE: NOvA's BpmAtTarget reads
//       double IY = BpmIntX[ind];  //yi[ind]
// i.e. the VERTICAL position is weighted with the HORIZONTAL intensity
// (the trailing comment shows yi was intended). We reproduce it so that
// posy agrees with NOvA spill by spill. Its size on real data is
// negligible: on the 2,735 spills of 2024-07-12 00:01-01:01 CDT the two
// weightings differ by a median of 2.4e-5 mm (max 2.6e-4 mm), against a
// 4 mm wide posy window.
//
// NOvA passes intensities already rescaled to POT units (see
// batchIntensitiesInPOT); that rescaling is a common factor that cancels
// in the weighted mean, so raw intensities give the identical result
// while keeping the position defined even when the POT reading is
// missing or zero (where NOvA's rescaling degenerates). Such spills fail
// the POT cut either way.
bool bpmAtTarget(const std::vector<double>& xp, const std::vector<double>& yp,
                 const std::vector<double>& xi, const std::vector<double>& yi,
                 double& posx, double& posy) {
    double IXtot = 0.0, IYtot = 0.0;
    double ix = 0.0, iy = 0.0;
    for (size_t ind = 0; ind < xi.size(); ++ind) {
        if (xi[ind] == 0.0 || yi[ind] == 0.0) continue;
        const double X = xp[ind];
        const double Y = yp[ind];
        double IX = xi[ind];
        double IY = xi[ind];   // NOvA quirk: horizontal intensity for the y weight (see above)
        IX = IX > 0 ? IX : 0;
        IY = IY > 0 ? IY : 0;
        IXtot += IX;  IYtot += IY;   // "sum intensities"
        ix += IX * X; iy += IY * Y;  // "intensity-weighted position value"
    }
    if (IXtot <= 0.0 || IYtot <= 0.0) return false;
    posx = ix / IXtot;   // "intensity-weighted average positions"
    posy = iy / IYtot;
    return true;
}

// NOvA's per-batch intensities in POT units ("convert batch by batch
// intensities in POT units"): each batch's BPM intensity times
// spillpot / sum(intensities), so the batches sum to the spill POT.
// NOvA's values are in units of 1e12 ("units will be e12"); here they
// are absolute protons, like spillpot. Left empty when spillpot or the
// intensity sum is not positive (undefined; NOvA would produce 0 or
// non-finite entries there).
std::vector<double> batchIntensitiesInPOT(const std::vector<double>& intensity, double spillpot) {
    double tot = 0.0;
    for (double v : intensity) tot += v;   // NOvA sums all entries, negatives included
    if (!(spillpot > 0.0) || !(tot > 0.0)) return {};
    std::vector<double> out;
    out.reserve(intensity.size());
    for (double v : intensity) out.push_back(v * spillpot / tot);
    return out;
}

// Port of GetStats [2]: "Calculate mean and standard deviation for
// points in TGraph to use as starting parameters values in Fit". The
// profile has already been inverted, so it is a Gaussian with a negative
// offset; "To get good starting param values, shift all profile values
// up by the most negative of all profile values". As in NOvA the shift
// is applied to the graph itself, so the subsequent fit runs on the
// shifted profile. Returns the total (shifted) charge.
double profileStats(TGraph& prof, double& mean, double& stdev) {
    mean = 0.0; stdev = 0.0;
    double minprof = 1e6;
    for (int i = 0; i < prof.GetN(); ++i) minprof = std::min(minprof, prof.GetY()[i]);
    double qtot = 0, qx = 0;
    for (int i = 0; i < prof.GetN(); ++i) {
        prof.SetPoint(i, prof.GetX()[i], prof.GetY()[i] - minprof);
        qtot += prof.GetY()[i];                     // "Total charge on all wires"
        qx += prof.GetY()[i] * prof.GetX()[i];      // "Total charge-weighted position value"
    }
    if (qtot == 0.0) return qtot;                   // "If total charge is 0, mean and stdev are 0"
    mean = qx / qtot;
    double var = 0;                                 // "Charge-weighted variance"
    for (int i = 0; i < prof.GetN(); ++i) {
        double d = prof.GetX()[i] - mean;
        var += prof.GetY()[i] * d * d;
    }
    var /= qtot;
    if (var > 0.0) stdev = std::sqrt(var);
    return qtot;
}

// Port of ProfileProjection + GetGaussFit [2]: Gaussian fit to one
// plane's 48 multiwire channels. IFBeam returns the wire readings in
// volts; NOvA divides them by Munits::millivolt (= 1e-3 V, since NOvA's
// Munits.h defines volt = kg m^2 / (A s^3) = 1 in its base units), i.e.
// works in mV — hence the x1000 below, which keeps NOvA's 100 mV charge
// threshold and its fit start values identical.
//   - wires at (ch - 24.5) * 0.5 mm: "Origin is halfway between channels
//     24 and 25";
//   - "Invert voltage to make it positive in Gaussian peak - may lead to
//     better fitting";
//   - dead channels dropped: "I WILL REMOVE DEAD CHANNELS (IE THOSE WITH
//     ZERO READING)"; no clamping to >= 0 ("INSTEAD FIT TO AN OFFSET").
// Returns the fitted sigma [mm] ("Standard deviation values are not
// extrapolated to target"), or NaN when NOvA's preconditions or
// physicality checks fail (NOvA then keeps its -9 m default).
double profileGaussWidth(const std::vector<double>& wires) {
    const double bad = std::numeric_limits<double>::quiet_NaN();
    if (wires.size() < 48) return bad;
    TGraph prof;
    for (int ch = 1; ch <= 48; ++ch) {
        double xpos = (ch - 24.5) * kWireSpacing;   // mm
        double qmv = wires[ch - 1] * 1000.0;        // V -> mV (NOvA's /Munits::millivolt)
        prof.SetPoint(prof.GetN(), xpos, -qmv);     // invert
    }
    for (int i = 0; i < prof.GetN(); ++i)
        if (prof.GetY()[i] == 0.0) { prof.RemovePoint(i); --i; }

    double statMean = 0, statStdev = 0;
    double area = profileStats(prof, statMean, statStdev);   // shifts profile

    // "Check we have enough points and enough charge for a reasonable fit":
    // > 5 live wires, > 100 mV total, stats mean within 20 wires (10 mm),
    // stats stdev within 12 wires (6 mm).
    if (!(prof.GetN() > 5 && std::abs(area) > 100.0 &&
          std::abs(statMean) < 20.0 * kWireSpacing &&
          std::abs(statStdev) < 12.0 * kWireSpacing))
        return bad;

    // "For some reason this is faster than if I use a TMath::Gaus"
    TF1 gausF("gausF",
        "[3]+([0]*exp((-1*(x-[1])*(x-[1]))/(2*[2]*[2]))/(sqrt(2*3.142)*[2]))");
    // Start values as in NOvA. Its offset start is "the avg value of the
    // negative readings", falling back to -3.0 ("a guess based on looking
    // at data"); because GetStats has already shifted the profile to be
    // non-negative, no negative readings remain and NOvA always uses -3.0.
    gausF.SetParameters(area, statMean, statStdev / 2.0, -3.0);
    if (prof.Fit(&gausF, "q", "", -24 * kWireSpacing, 24 * kWireSpacing) != 0)
        return bad;

    // "Are fit results physical?"
    if (gausF.GetParameter(0) < 10.0) return bad;
    if (std::abs(gausF.GetParameter(1)) > 24 * kWireSpacing) return bad;
    if (gausF.GetParameter(2) > 24 * kWireSpacing || gausF.GetParameter(2) < 0.0) return bad;
    return gausF.GetParameter(2);   // sigma [mm]
}

// Calibrated total horn current from the four stripline readbacks
// (E:NSLINA-D, "current on NuMI strip line A..D", each ~1/4 of the total),
// verbatim from IFDBSpillInfo::produce [2]: "Updated numbers from Jim
// Hylen, email Oct 4, 2013".
double calibratedHornCurrent(double nslina, double nslinb, double nslinc, double nslind) {
    double ihorn = 0.;
    ihorn += (nslina - (+0.01)) / 0.9951;
    ihorn += (nslinb - (-0.14)) / 0.9957;
    ihorn += (nslinc - (-0.05)) / 0.9965;
    ihorn += (nslind - (-0.07)) / 0.9945;
    return ihorn;
}

// ---------------------------------------------------------------------
// Main entry point (ROOT macro convention: function name == file name).
//
//   root -l -b -q 'get_data.cpp("2024-07-12T00:01:00-05:00","2024-07-12T01:01:00-05:00")'
//   root -l -b -q 'get_data.cpp("<t0>","<t1>","P11-12")'   # older NOvA epoch
//
// Times are ISO-8601 with an explicit UTC offset (-05:00 = CDT).
// novaEpoch selects NOvA's period-dependent position-y/width cuts
// (see novaCutsForEpoch); the default "P13+" covers all 2x2 running.
// ---------------------------------------------------------------------
int get_data(const std::string& min_time_iso = "2024-07-12T00:01:00-05:00",
             const std::string& max_time_iso = "2024-07-12T01:01:00-05:00",
             const std::string& novaEpoch = "P13+") {

    const NOvABeamCuts cuts = novaCutsForEpoch(novaEpoch);

    // -----------------------------------------------------------------
    // The ACNET devices to pull from the IFBeam archive. All are real
    // ACNET device names except E:NSLIN, which is a pseudo-device
    // computed below. Per-device notes follow [1] Sec. 4 and the inline
    // documentation of [2].
    // -----------------------------------------------------------------
    DeviceMap deviceMap = {
        //_________POT (protons on target, per spill)_________
        // Two toroidal beam-current transformers encircle the primary
        // proton beam [1]: TOR101 sits in the upstream NuMI beamline,
        // TORTGT ~10 m before the target. The trailing "D" selects the
        // $A9-datalogged readback of each toroid; [2] describes
        // "E:TRTGTD -- 16-bit integer of intensity (POT) on NuMI target
        // toroid" and "E:TR101D -- 16-bit integer of intensity (POT) on
        // NuMI toroid 101". IFBeam units are "E12" (readings in units of
        // 1e12 protons; unitToFactor() applies the factor, so values
        // here are absolute protons/spill, ~5e13 at nominal intensity).
        // [2] uses TRTGTD as the primary POT device and falls back to
        // TR101D when TRTGTD is missing or unphysically small.
        {"E:TRTGTD", {}},   // POT at target toroid TORTGT (primary POT device)
        {"E:TOR101", {}},   // upstream toroid TOR101, direct readback (not used by NOvA)
        {"E:TR101D", {}},   // upstream toroid TOR101, datalogged readback (NOvA POT backup)

        //_________Horn current_______________________________
        // The NuMI focusing-horn current (~200 kA nominal [1]) is
        // delivered through parallel stripline conductors; NSLINA-D are
        // the four stripline current transductor readbacks (each ~1/4
        // of the total, in kA, NEGATIVE in both horn polarities). [2]
        // documents them as "current on NuMI strip line A/B/C/D"; the
        // calibrated sum (calibratedHornCurrent) is the total.
        {"E:NSLINA", {}},   // stripline A current [kA]
        {"E:NSLINB", {}},   // stripline B current [kA]
        {"E:NSLINC", {}},   // stripline C current [kA]
        {"E:NSLIND", {}},   // stripline D current [kA]
        // Pseudo-device: NOT fetched from IFBeam. Filled below with the
        // calibrated total horn current computed from NSLINA-D.
        {"E:NSLIN", {}},    // total horn current [kA] (computed here)

        //________Horn polarity_______________________________
        // Horn direction readback. Validated against NOvA's epoch
        // definitions (see "Horn polarity" in the header): about
        // -750 V in FHC; about +742 V in RHC in 2017-2019 but only
        // +19 V in RHC since the 2023-24 shutdown. E:HRNDIR > 0 <=> RHC
        // throughout — the convention of NOvA's 2016 (later removed)
        // logic.
        {"E:HRNDIR", {}},

        //__________Beam position (BPMs)______________________
        // Capacitive beam position monitors in the primary line [1]:
        // station "121" is upstream in the NuMI line, "TGT" is the last
        // station before the target. H/V prefix = horizontal/vertical
        // plane; "P" = position [mm], "I" = the intensity seen by the
        // BPM (used to weight/validate the position measurements).
        // The "[]" suffix requests the full array: 7 elements per
        // spill. Per [2]: element 0 is "an average value of a few
        // batches (often 2,3,4) used for auto-tuning, so we should
        // disregard it"; elements 1-6 are the per-Booster-batch
        // readings (a NuMI spill is 6 Booster batches).
        {"E:HPTGT[]", {}},  // horizontal position at target BPM [mm], 7 elements
        {"E:HITGT[]", {}},  // intensity, horizontal plane, target BPM
        {"E:HP121[]", {}},  // horizontal position at station-121 BPM [mm]
        {"E:VPTGT[]", {}},  // vertical position at target BPM [mm]
        {"E:VITGT[]", {}},  // intensity, vertical plane, target BPM
        {"E:VP121[]", {}},  // vertical position at station-121 BPM [mm]

        //__________Beam width (profile monitor)______________
        // Multiwire profile monitor at the target ("MTGTDS") [1]. The
        // 216-element array readback contains, per [2]
        // (ProfileProjection): the 48 horizontal-wire voltages at
        // indices [103..150] and the 48 vertical-wire voltages at
        // [151..198] (remaining elements are header/housekeeping
        // words). A Gaussian fit to each plane's wire profile yields
        // the beam width (sigma) in mm at the target.
        {"E:MTGTDS[]", {}}
    };

    // -----------------------------------------------------------------
    // Fetch every device's archive for the requested window.
    // Failures (network, device missing from the archive) skip that
    // device and continue — the same per-category tolerance NOvA adopted
    // in 2024 so that one missing device does not blank the others.
    // -----------------------------------------------------------------
    for (auto& [deviceName, beamInfo] : deviceMap) {
        beamInfo.clear();
        if (deviceName == "E:NSLIN") continue;   // pseudo-device, computed below
        std::string url_device = createUrl(deviceName, min_time_iso, max_time_iso);
        std::cout << url_device << "\n";
        try {
            beamInfo = retrieveInfoFromDataBase(url_device);
        } catch (const std::exception& e) {
            std::cerr << "Error on device '" << deviceName << "': " << e.what() << "\n";
            // continue with the other devices instead of crashing everything
        }
        std::cout << deviceName << ": " << beamInfo.size() << " spills\n";
    }

    // -----------------------------------------------------------------
    // Total horn current from the four stripline readbacks.
    //
    // For every NSLINA record, the B, C and D records are looked up
    // within beamMatchDT; the calibrated sum is stored under the pseudo-
    // device E:NSLIN only when all four exist for that spill. (An earlier
    // version kept B/C/D in variables shared across spills, so a missing
    // record silently reused the previous spill's value, and it skipped
    // any exactly-zero reading — which NOvA does not do, and which would
    // hide genuine horn-off readings. Both are fixed here.)
    // -----------------------------------------------------------------
    for (const auto& [timeA, valuesA] : deviceMap["E:NSLINA"]) {
        const std::vector<double>* b = findMatch(deviceMap["E:NSLINB"], timeA, beamMatchDT);
        const std::vector<double>* c = findMatch(deviceMap["E:NSLINC"], timeA, beamMatchDT);
        const std::vector<double>* d = findMatch(deviceMap["E:NSLIND"], timeA, beamMatchDT);
        if (valuesA.empty() || !b || b->empty() || !c || c->empty() || !d || d->empty())
            continue;   // NOvA: the four-stripline query throws -> no horn current
        deviceMap["E:NSLIN"][timeA] = {calibratedHornCurrent(valuesA.front(), b->front(),
                                                             c->front(), d->front())};
    }

    // -----------------------------------------------------------------
    // Output 1: ROOT TTree, one entry per spill.
    // The device with the most spills defines the reference timeline;
    // every other device is matched to it within beamMatchDT. Unmatched
    // scalars are NaN, unmatched arrays are empty vectors.
    // -----------------------------------------------------------------
    std::string refDevice;
    size_t refCount = 0;
    for (const auto& [name, info] : deviceMap)
        if (info.size() > refCount) { refCount = info.size(); refDevice = name; }
    if (refCount == 0) {
        std::cerr << "No data retrieved for any device; no output written.\n";
        return 1;
    }
    std::cout << "Reference timeline: " << refDevice << " (" << refCount << " spills)\n";

    // Filenames carry the queried time range. ':' is illegal in filenames
    // on Windows/OneDrive, so it is stripped from the ISO timestamps:
    // "2024-07-12T00:01:00-05:00" -> "2024-07-12T000100-0500"
    auto sanitizeTime = [](std::string s) {
        s.erase(std::remove(s.begin(), s.end(), ':'), s.end());
        return s;
    };
    const std::string fileStem = "get_data_" + sanitizeTime(min_time_iso) +
                                 "_to_" + sanitizeTime(max_time_iso);
    const std::string rootFileName = fileStem + ".root";
    const std::string csvFileName  = fileStem + ".csv";

    // "E:HPTGT[]" -> branch/column name "HPTGT"
    auto branchName = [](std::string name) {
        if (name.rfind("E:", 0) == 0) name = name.substr(2);
        while (!name.empty() && (name.back() == '[' || name.back() == ']')) name.pop_back();
        return name;
    };

    // Devices whose spills only ever carry one value get a plain double
    // branch; the rest get a std::vector<double> branch.
    std::map<std::string, bool> isScalar;
    for (const auto& [name, info] : deviceMap) {
        size_t maxSize = 0;
        for (const auto& [t, v] : info) maxSize = std::max(maxSize, v.size());
        isScalar[name] = (maxSize <= 1);
    }

    TFile fout(rootFileName.c_str(), "RECREATE");
    TTree tree("beam", "NuMI IFBeam per-spill data");
    double spillTime = 0.0;
    std::string spillISO;
    tree.Branch("time", &spillTime);
    tree.Branch("time_iso", &spillISO);
    std::map<std::string, double> scalarBuf;
    std::map<std::string, std::vector<double>> vectorBuf;
    for (const auto& [name, info] : deviceMap) {
        const std::string bname = branchName(name);
        if (isScalar[name]) tree.Branch(bname.c_str(), &scalarBuf[name]);
        else                tree.Branch(bname.c_str(), &vectorBuf[name]);
    }

    // -----------------------------------------------------------------
    // Parallel NOvA-style "quality" tree, entry-aligned with "beam" so
    // it can be attached as a friend tree. Evaluates goodbeam cuts 1-5
    // (see the "Beam-quality cuts" header comment; cut 6 needs a
    // detector trigger stream and is not evaluated). Nothing is
    // filtered — failing spills are written with goodbeam = 0.
    // -----------------------------------------------------------------
    TTree qtree("quality", "NOvA-style per-spill goodbeam evaluation (cuts 1-5)");
    const double dNaN = std::numeric_limits<double>::quiet_NaN();
    double q_spilltime, q_pot, q_horn, q_horndir, q_posx, q_posy, q_wx, q_wy;
    std::vector<double> q_bposx, q_bposy, q_intx, q_inty;
    Bool_t p_pot, p_horn, p_posx, p_posy, p_wx, p_wy, q_good, q_is0hc, q_isrhc;
    qtree.Branch("time", &spillTime);         // reference-timeline time [s] (as "beam")
    qtree.Branch("time_iso", &spillISO);
    qtree.Branch("spilltime", &q_spilltime);  // NOvA's "ifdbtime": TRTGTD record time, else TR101D [s]
    qtree.Branch("spillpot", &q_pot);         // NOvA POT logic (TRTGTD w/ TR101D fallback) [protons]
    qtree.Branch("hornI", &q_horn);           // calibrated NSLIN sum [kA]
    qtree.Branch("hornDir", &q_horndir);      // E:HRNDIR readback [V]
    qtree.Branch("isRHC", &q_isrhc);          // E:HRNDIR > 0 (validated convention)
    qtree.Branch("is0HC", &q_is0hc);          // |hornI| < 1 kA (horn off)
    qtree.Branch("posx", &q_posx);            // intensity-weighted, extrapolated to target [mm]
    qtree.Branch("posy", &q_posy);
    qtree.Branch("widthx", &q_wx);            // multiwire Gaussian sigma [mm]
    qtree.Branch("widthy", &q_wy);
    qtree.Branch("bposx", &q_bposx);          // per-batch target positions [mm] (NOvA bposx; 0 = skipped batch)
    qtree.Branch("bposy", &q_bposy);
    qtree.Branch("intx", &q_intx);            // per-batch intensities in POT [protons] (NOvA intx, x1e12)
    qtree.Branch("inty", &q_inty);
    qtree.Branch("pass_pot", &p_pot);
    qtree.Branch("pass_horn", &p_horn);
    qtree.Branch("pass_posx", &p_posx);
    qtree.Branch("pass_posy", &p_posy);
    qtree.Branch("pass_widthx", &p_wx);
    qtree.Branch("pass_widthy", &p_wy);
    qtree.Branch("goodbeam", &q_good);

    std::ofstream qcsv(fileStem + "_quality.csv");
    // The first 15 columns keep their original order so older readers
    // keep working; the reconciliation added the last four.
    qcsv << "time_s,time_iso,spillpot,hornI,posx,posy,widthx,widthy,"
            "pass_pot,pass_horn,pass_posx,pass_posy,pass_widthx,pass_widthy,goodbeam,"
            "spilltime_s,hornDir,isRHC,is0HC\n";
    qcsv << std::setprecision(15);
    long long nGood = 0, nRHC = 0, n0HC = 0;
    long long nFail[6] = {0, 0, 0, 0, 0, 0};

    for (const auto& [t, refValues] : deviceMap[refDevice]) {
        spillTime = t;
        spillISO = toISO8601(t);
        for (const auto& [name, info] : deviceMap) {
            const std::vector<double>* match = findMatch(info, t, beamMatchDT);
            if (isScalar[name]) {
                scalarBuf[name] = (match && !match->empty())
                                      ? match->front()
                                      : std::numeric_limits<double>::quiet_NaN();
            } else {
                vectorBuf[name] = match ? *match : std::vector<double>{};
            }
        }
        tree.Fill();

        // ---- NOvA-style quality evaluation for this spill ----

        // 1. POT. NOvA queries TRTGTD and TR101D separately ("This
        //    ensures we can use TRTGTD still if TR101D doesn't exist"):
        //      "If we found TRTGTD in IFDB": POT = TRTGTD, replaced by
        //        TR101D when TR101D is nonzero and TRTGTD < 0.02e12;
        //      "only reach here if TRTGTD is missing": POT = TR101D.
        //    The record time of whichever toroid was found first is
        //    NOvA's "ifdbtime" — the spill time the timing cut compares
        //    with the trigger. Negative readings are clamped to 0.
        auto itTRTGTD = findMatchIt(deviceMap["E:TRTGTD"], t, beamMatchDT);
        auto itTR101D = findMatchIt(deviceMap["E:TR101D"], t, beamMatchDT);
        const bool haveTRTGTD = itTRTGTD != deviceMap["E:TRTGTD"].end() && !itTRTGTD->second.empty();
        const bool haveTR101D = itTR101D != deviceMap["E:TR101D"].end() && !itTR101D->second.empty();
        q_pot = dNaN;
        q_spilltime = dNaN;
        if (haveTRTGTD) {
            q_spilltime = itTRTGTD->first;
            q_pot = itTRTGTD->second.front();
            if (haveTR101D && itTR101D->second.front() != 0.0 && q_pot < 0.02e12)
                q_pot = itTR101D->second.front();
        } else if (haveTR101D) {
            q_spilltime = itTR101D->first;
            q_pot = itTR101D->second.front();
        }
        if (q_pot < 0.0) q_pot = 0.0;   // NaN compares false and stays NaN
        p_pot = !std::isnan(q_pot) && q_pot > cuts.minPOT;

        // 2. Horn current and polarity.
        q_horn = scalarBuf["E:NSLIN"];
        p_horn = q_horn > cuts.minHornI && q_horn < cuts.maxHornI;
        q_is0hc = std::fabs(q_horn) < cuts.hornI0HC;   // "Horn Off, flagged when abs(ihorn) < 1"
        q_horndir = scalarBuf["E:HRNDIR"];
        q_isrhc = q_horndir > 0.0;                     // NaN (missing) -> false

        // 3-4. Position at target: per-batch BPM extrapolation and the
        //      intensity-weighted mean (NOvA BpmProjection + BpmAtTarget).
        std::vector<double> xp, yp, xi, yi;
        bpmProjection(vectorBuf["E:VP121[]"], vectorBuf["E:HP121[]"],
                      vectorBuf["E:VPTGT[]"], vectorBuf["E:HPTGT[]"],
                      vectorBuf["E:VITGT[]"], vectorBuf["E:HITGT[]"],
                      xp, yp, xi, yi);
        q_bposx = xp;
        q_bposy = yp;
        q_intx = batchIntensitiesInPOT(xi, q_pot);
        q_inty = batchIntensitiesInPOT(yi, q_pot);
        q_posx = dNaN; q_posy = dNaN;
        bpmAtTarget(xp, yp, xi, yi, q_posx, q_posy);
        p_posx = q_posx > cuts.minPosX && q_posx < cuts.maxPosX;
        p_posy = q_posy > cuts.minPosY && q_posy < cuts.maxPosY;

        // 5. Width: Gaussian sigma of each multiwire plane
        //    (MTGTDS elements 103-150 horizontal, 151-198 vertical;
        //    NOvA's ProfileProjection indexes PMTGT[103+ind] and
        //    PMTGT[103+48+ind] for ind = 0..47).
        const std::vector<double>& mw = vectorBuf["E:MTGTDS[]"];
        q_wx = dNaN; q_wy = dNaN;
        if (mw.size() >= 199) {
            q_wx = profileGaussWidth({mw.begin() + 103, mw.begin() + 151});
            q_wy = profileGaussWidth({mw.begin() + 151, mw.begin() + 199});
        }
        p_wx = q_wx > cuts.minWidthX && q_wx < cuts.maxWidthX;
        p_wy = q_wy > cuts.minWidthY && q_wy < cuts.maxWidthY;

        // "Apply good beam cuts" (cuts 1-5; cut 6 needs a trigger time).
        q_good = p_pot && p_horn && p_posx && p_posy && p_wx && p_wy;
        if (q_good) ++nGood;
        if (q_isrhc) ++nRHC;
        if (q_is0hc) ++n0HC;
        nFail[0] += !p_pot; nFail[1] += !p_horn; nFail[2] += !p_posx;
        nFail[3] += !p_posy; nFail[4] += !p_wx; nFail[5] += !p_wy;
        qtree.Fill();
        qcsv << spillTime << ',' << spillISO << ',' << q_pot << ',' << q_horn << ','
             << q_posx << ',' << q_posy << ',' << q_wx << ',' << q_wy << ','
             << int(p_pot) << ',' << int(p_horn) << ',' << int(p_posx) << ','
             << int(p_posy) << ',' << int(p_wx) << ',' << int(p_wy) << ','
             << int(q_good) << ',' << q_spilltime << ',' << q_horndir << ','
             << int(q_isrhc) << ',' << int(q_is0hc) << '\n';
    }
    tree.Write();
    qtree.Write();
    qcsv.close();

    // -----------------------------------------------------------------
    // Output 3: the cut values actually applied, so plots and reports
    // never drift from the evaluation. Stored as TNamed/TParameter
    // objects in the ROOT file and as a JSON sidecar.
    // -----------------------------------------------------------------
    const std::vector<std::pair<std::string, double>> cutValues = {
        {"minPOT", cuts.minPOT},       {"minHornI", cuts.minHornI},
        {"maxHornI", cuts.maxHornI},   {"hornI0HC", cuts.hornI0HC},
        {"minPosX", cuts.minPosX},     {"maxPosX", cuts.maxPosX},
        {"minPosY", cuts.minPosY},     {"maxPosY", cuts.maxPosY},
        {"minWidthX", cuts.minWidthX}, {"maxWidthX", cuts.maxWidthX},
        {"minWidthY", cuts.minWidthY}, {"maxWidthY", cuts.maxWidthY},
        {"maxDeltaT", cuts.maxDeltaT}};
    TNamed("cut_epoch", cuts.epoch.c_str()).Write();
    for (const auto& [name, value] : cutValues)
        TParameter<double>(("cut_" + name).c_str(), value).Write();
    json cutJson;
    cutJson["epoch"] = cuts.epoch;
    cutJson["source"] = "NOvA IFDBSpillInfo.fcl standard_ifdbspillinfo (novasoft main), "
                        "converted to absolute protons, kA, mm, s";
    for (const auto& [name, value] : cutValues) cutJson[name] = value;
    std::ofstream(fileStem + "_cuts.json") << cutJson.dump(2) << '\n';

    std::cout << "Wrote " << tree.GetEntries() << " spills to " << rootFileName
              << " (TTrees 'beam' + 'quality', cut values as cut_* objects)\n";
    std::cout << "NOvA epoch " << cuts.epoch << ": goodbeam (cuts 1-5) " << nGood << " / "
              << qtree.GetEntries() << " spills pass; failures by cut: POT " << nFail[0]
              << ", horn " << nFail[1] << ", posx " << nFail[2] << ", posy " << nFail[3]
              << ", widthx " << nFail[4] << ", widthy " << nFail[5] << "\n";
    std::cout << "Polarity: " << nRHC << " spills with E:HRNDIR > 0 (RHC), "
              << n0HC << " horn-off (|hornI| < 1 kA)\n";
    std::cout << "Quality CSV: " << fileStem << "_quality.csv; cuts: " << fileStem << "_cuts.json\n";
    fout.Close();

    // -----------------------------------------------------------------
    // Output 2: CSV in long/tidy format, one row per (device, spill,
    // array index) — no time matching, exactly what came back per device.
    // -----------------------------------------------------------------
    std::ofstream csv(csvFileName);
    csv << "device,time_s,time_iso,index,value\n";
    csv << std::setprecision(15);
    long long nRows = 0;
    for (const auto& [name, info] : deviceMap) {
        for (const auto& [t, values] : info) {
            const std::string iso = toISO8601(t);
            for (size_t i = 0; i < values.size(); ++i, ++nRows)
                csv << name << ',' << t << ',' << iso << ',' << i << ',' << values[i] << '\n';
        }
    }
    csv.close();
    std::cout << "Wrote " << nRows << " rows to " << csvFileName << "\n";

    return 0;
}
